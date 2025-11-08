#!/usr/bin/env python3
"""
Per-variant Fisher / Barnard / Bayesian pipeline (uses provided gnomAD_Homozygotes if present)
Fixed Barnard implementation for large samples
"""
import os
import gc
import re
import math
import warnings
from functools import partial
from concurrent.futures import ThreadPoolExecutor
from threading import Lock

import pandas as pd
import numpy as np
from scipy import stats
from statsmodels.stats.multitest import multipletests
from tqdm import tqdm

warnings.filterwarnings('ignore')

# Barnard exact test availability check
try:
    from scipy.stats import barnard_exact
    _has_barnard = True
except Exception:
    _has_barnard = False

# ----------------------------
# CONFIG - edit paths / params here
# ----------------------------
CONFIG = {
    'input_file': 'variant_analysis_full.csv',            
    'output_variant_fisher': 'variant_fisher.csv',
    'output_variant_barnard': 'variant_barnard.csv',
    'output_variant_bayesian': 'variant_bayesian.csv',
    'output_skipped': 'skipped_variants.csv',
    'total_cases': 62,
    # Bayesian sampling - keep moderate for performance; increase if you need more precision
    'bayesian_nsam': 1000000,
    'bayesian_prior_alpha': 0.5,
    'bayesian_prior_beta': 0.5,
    'significance_level': 0.05,
    'gnomad_ac_col': 'gnomAD_Allele_Count',
    'gnomad_an_col': 'gnomAD_Allele_Number',
    'gnomad_hom_col': 'gnomAD_Homozygotes',   
    # Parallel / chunking
    'chunk_size': 12,   
    'n_workers': 12,    
    'flush_every_chunk': True,
    # Barnard test limits
    'barnard_max_total': 10000,  # Skip Barnard if total sample size exceeds this
    'barnard_downsample_to': 5000,  # Downsample controls to this size for Barnard
    'barnard_timeout': 180  # Maximum seconds for Barnard test
}

# ----------------------------
# Utilities
# ----------------------------
re_split_pattern = re.compile(r'[;,|\s]+')


def parse_sample_count(files_field):
    """Return number of sample IDs in SRR_Files-like field (fast)."""
    if pd.isna(files_field) or str(files_field).strip() == '':
        return 0
    s = str(files_field)
    parts = [p.strip() for p in re_split_pattern.split(s) if p.strip()]
    return len(parts)


def calculate_odds_ratio_ci(a, b, c, d, alpha=0.05):
    """Compute OR and Woolf CI with Haldane-Anscombe correction."""
    a_adj = a + 0.5
    b_adj = b + 0.5
    c_adj = c + 0.5
    d_adj = d + 0.5
    or_val = (a_adj * d_adj) / (b_adj * c_adj)
    se_log_or = math.sqrt(1.0/a_adj + 1.0/b_adj + 1.0/c_adj + 1.0/d_adj)
    z = stats.norm.ppf(1 - alpha/2)
    log_or = math.log(or_val)
    ci_lower = math.exp(log_or - z * se_log_or)
    ci_upper = math.exp(log_or + z * se_log_or)
    return or_val, ci_lower, ci_upper


def expected_control_carriers_from_ac_an(ac, an):
    """
    Convert AC/AN -> expected number of control individuals carrying variant.
    Uses Hardy-Weinberg conversion: p = AC/AN; N = AN/2; expected = N*(2*p - p^2)
    Return np.nan if inputs invalid.
    """
    if pd.isna(ac) or pd.isna(an):
        return np.nan
    try:
        acf = float(ac)
        anf = float(an)
    except Exception:
        return np.nan
    if anf <= 0:
        return np.nan
    p = acf / anf
    n_indiv = anf / 2.0
    expected = n_indiv * (2.0 * p - p * p)
    return float(expected)


def compute_barnard_with_downsampling(a, b, c, d, config):
    """
    Compute Barnard test with intelligent downsampling for large control groups.
    Returns p-value and method used ('exact', 'downsampled', 'chi2_approx', or 'failed')
    """
    total_n = a + b + c + d
    ctrl_n = c + d
    
    # If total sample size is manageable, use exact test
    if total_n <= config.get('barnard_max_total', 10000):
        try:
            import signal
            from contextlib import contextmanager
            
            # Timeout handler
            @contextmanager
            def timeout(seconds):
                def timeout_handler(signum, frame):
                    raise TimeoutError()
                old_handler = signal.signal(signal.SIGALRM, timeout_handler)
                signal.alarm(seconds)
                try:
                    yield
                finally:
                    signal.alarm(0)
                    signal.signal(signal.SIGALRM, old_handler)
            
            # Try exact test with timeout (Unix-like systems only)
            try:
                with timeout(config.get('barnard_timeout', 10)):
                    table = [[a, b], [c, d]]
                    res = barnard_exact(table, alternative='two-sided')
                    return res.pvalue, 'exact'
            except:
                # Fallback for Windows or timeout
                table = [[a, b], [c, d]]
                res = barnard_exact(table, alternative='two-sided')
                return res.pvalue, 'exact'
                
        except Exception:
            pass
    
    # If control group is large, use downsampling
    downsample_size = config.get('barnard_downsample_to', 5000)
    if ctrl_n > downsample_size:
        # Maintain the same carrier proportion when downsampling
        carrier_prop = c / ctrl_n if ctrl_n > 0 else 0
        c_down = int(round(carrier_prop * downsample_size))
        d_down = downsample_size - c_down
        
        try:
            table_down = [[a, b], [c_down, d_down]]
            res = barnard_exact(table_down, alternative='two-sided')
            return res.pvalue, 'downsampled'
        except Exception:
            pass
    
    # Fallback to chi-squared approximation for very large samples
    try:
        table = [[a, b], [c, d]]
        # Use chi2 with Yates correction for 2x2 table
        chi2, p_chi2, _, _ = stats.chi2_contingency(table, correction=True)
        return p_chi2, 'chi2_approx'
    except Exception:
        return np.nan, 'failed'


# ----------------------------
# Thread-safe file append helper
# ----------------------------
file_locks = {
    'fisher': Lock(),
    'barnard': Lock(),
    'bayes': Lock(),
    'skipped': Lock()
}


def append_results_to_csv(df, filepath, write_header, lock):
    """Thread-safe append DataFrame to CSV file"""
    if df is None or df.empty:
        return
    # ensure directory exists
    d = os.path.dirname(filepath)
    if d and not os.path.exists(d):
        os.makedirs(d, exist_ok=True)
    with lock:
        mode = 'w' if write_header else 'a'
        df.to_csv(filepath, mode=mode, header=write_header, index=False)


def initialize_output_files(config):
    """Remove existing outputs at start (called before multithreading)."""
    for filepath in [config['output_variant_fisher'],
                     config['output_variant_barnard'],
                     config['output_variant_bayesian'],
                     config['output_skipped']]:
        try:
            if os.path.exists(filepath):
                os.remove(filepath)
        except Exception:
            pass


# ----------------------------
# Per-variant processing (thread-safe)
# ----------------------------
def process_single_variant(row_data, config):
    """
    row_data: tuple (idx, variant_id, files_field, file_count, ac_val, an_val, hom_val)
    """
    idx, variant_id, files_field, file_count, ac_val, an_val, hom_val = row_data

    # Cases carriers
    try:
        carriers = parse_sample_count(files_field)
    except Exception:
        carriers = 0
    if carriers == 0:
        try:
            carriers = int(float(file_count or 0))
        except Exception:
            carriers = 0
    a = int(max(carriers, 0))
    b = int(max(config['total_cases'] - a, 0))

    # Strict requirement: numeric AC & AN and AN > 0
    if pd.isna(ac_val) or pd.isna(an_val):
        return {'status': 'skipped', 'skipped': {
            'Variant': variant_id, 'idx': idx, 'reason': 'missing_AC_or_AN',
            'gnomAD_AC_raw': ac_val, 'gnomAD_AN_raw': an_val, 'Cases_Carriers': a
        }}
    try:
        acf = float(ac_val)
        anf = float(an_val)
    except Exception:
        return {'status': 'skipped', 'skipped': {
            'Variant': variant_id, 'idx': idx, 'reason': 'AC_AN_not_numeric',
            'gnomAD_AC_raw': ac_val, 'gnomAD_AN_raw': an_val, 'Cases_Carriers': a
        }}
    if anf <= 0:
        return {'status': 'skipped', 'skipped': {
            'Variant': variant_id, 'idx': idx, 'reason': 'AN_leq_zero',
            'gnomAD_AC_raw': ac_val, 'gnomAD_AN_raw': an_val, 'Cases_Carriers': a
        }}

    # Control sample size (individuals)
    ctrl_N_indiv = int(round(anf / 2.0))
    if ctrl_N_indiv <= 0:
        return {'status': 'skipped', 'skipped': {
            'Variant': variant_id, 'idx': idx, 'reason': 'Control_N_zero',
            'gnomAD_AC_raw': acf, 'gnomAD_AN_raw': anf, 'Cases_Carriers': a
        }}

    # If homozygote count is provided and numeric, use it to estimate control carriers
    c = None
    if hom_val is not None and not pd.isna(hom_val):
        try:
            homf = float(hom_val)
            hom_int = int(round(max(homf, 0)))
            # cap homozygotes to number of individuals
            if hom_int > ctrl_N_indiv:
                hom_int = ctrl_N_indiv
            # estimate heterozygotes from allele counts: het_est = AC - 2*homo
            het_est = acf - 2.0 * hom_int
            if het_est < 0:
                het_est = 0.0
            # each heterozygote contributes 1 alt allele
            het_int = int(round(het_est))
            controls_est_carriers = hom_int + het_int
            # cap total carriers to number of individuals
            if controls_est_carriers > ctrl_N_indiv:
                controls_est_carriers = ctrl_N_indiv
            c = int(round(controls_est_carriers))
        except Exception:
            c = None

    # fallback: if we couldn't compute c from homozygotes, use expected HWE-based method
    if c is None:
        expected_ctrl = expected_control_carriers_from_ac_an(acf, anf)
        if pd.isna(expected_ctrl):
            return {'status': 'skipped', 'skipped': {
                'Variant': variant_id, 'idx': idx, 'reason': 'expected_ctrl_nan',
                'gnomAD_AC_raw': acf, 'gnomAD_AN_raw': anf, 'Cases_Carriers': a
            }}
        try:
            c = int(round(float(expected_ctrl)))
        except Exception:
            return {'status': 'skipped', 'skipped': {
                'Variant': variant_id, 'idx': idx, 'reason': 'expected_ctrl_not_int',
                'gnomAD_AC_raw': acf, 'gnomAD_AN_raw': anf, 'Cases_Carriers': a
            }}

    # ensure bounds
    if c < 0:
        c = 0
    if c > ctrl_N_indiv:
        c = ctrl_N_indiv
    d = int(max(ctrl_N_indiv - c, 0))

    # compute OR and CI
    try:
        orv, cil, cir = calculate_odds_ratio_ci(a, b, c, d)
    except Exception:
        orv, cil, cir = np.nan, np.nan, np.nan

    base = {
        'Variant': variant_id,
        'Cases_Carriers': a,
        'Cases_NonCarriers': b,
        'Controls_Est_Carriers': c,
        'Controls_NonCarriers': d,
        'gnomAD_AC': acf,
        'gnomAD_AN': anf,
        'gnomAD_Homozygotes': hom_val,
        'Control_N_indiv': ctrl_N_indiv,
        'Odds_Ratio': orv,
        'CI_Lower': cil,
        'CI_Upper': cir
    }

    # Fisher
    try:
        table = [[a, b], [c, d]]
        _, p_fisher = stats.fisher_exact(table, alternative='two-sided')
    except Exception:
        p_fisher = np.nan
    fisher_res = {**base, 'P_Value': p_fisher}

    # Barnard with intelligent handling for large samples
    p_barnard = np.nan
    barnard_method = 'not_available'
    if _has_barnard:
        p_barnard, barnard_method = compute_barnard_with_downsampling(a, b, c, d, config)
    
    barnard_res = {
        **base, 
        'P_Value': p_barnard, 
        'Barnard_available': _has_barnard,
        'Barnard_method': barnard_method
    }

    # Bayesian (Beta posteriors on carrier probability)
    nsam = int(config.get('bayesian_nsam', 1000000))
    alpha_prior = config.get('bayesian_prior_alpha', 0.5)
    beta_prior = config.get('bayesian_prior_beta', 0.5)
    # deterministic but different RNG per variant index
    rng = np.random.default_rng(42 + (idx % 100000))

    try:
        p_cases = rng.beta(a + alpha_prior, b + beta_prior, nsam)
        p_ctrls = rng.beta(c + alpha_prior, d + beta_prior, nsam)
        p_cases = np.clip(p_cases, 1e-10, 1 - 1e-10)
        p_ctrls = np.clip(p_ctrls, 1e-10, 1 - 1e-10)
        odds_cases = p_cases / (1 - p_cases)
        odds_ctrls = p_ctrls / (1 - p_ctrls)
        or_samples = odds_cases / odds_ctrls
        posterior_median_or = float(np.median(or_samples))
        or_ci_l = float(np.percentile(or_samples, 2.5))
        or_ci_u = float(np.percentile(or_samples, 97.5))
        posterior_prob_gt1 = float(np.mean(or_samples > 1.0))
        bayes_p = float(2.0 * min(posterior_prob_gt1, 1.0 - posterior_prob_gt1))
        # free memory
        del p_cases, p_ctrls, odds_cases, odds_ctrls, or_samples
    except Exception:
        posterior_median_or = np.nan
        or_ci_l = np.nan
        or_ci_u = np.nan
        posterior_prob_gt1 = np.nan
        bayes_p = np.nan

    bayes_res = {
        **base,
        'Posterior_Median_OR': posterior_median_or,
        'OR_CI_Lower': or_ci_l,
        'OR_CI_Upper': or_ci_u,
        'Posterior_Prob_OR_gt_1': posterior_prob_gt1,
        'Bayes_p': bayes_p
    }

    return {'status': 'processed', 'fisher': fisher_res, 'barnard': barnard_res, 'bayes': bayes_res}


# ----------------------------
# Prepare minimal variant data (reduces memory)
# ----------------------------
def prepare_variant_data(df, config):
    """Return list of tuples: (idx, variant_id, files_field, file_count, ac_val, an_val, hom_val)"""
    ac_col = config['gnomad_ac_col']
    an_col = config['gnomad_an_col']
    hom_col = config.get('gnomad_hom_col', 'gnomAD_Homozygotes')
    prepared = []
    for idx, r in df.iterrows():
        variant_id = r.get('Variant') or r.get('variant') or f"row{idx}"
        files_field = r.get('SRR_Files', '')
        file_count = r.get('File_Count', 0)
        # robust numeric parsing
        ac_raw = r.get(ac_col, None)
        an_raw = r.get(an_col, None)
        hom_raw = r.get(hom_col, None) if hom_col in r else None
        ac_val = pd.to_numeric(ac_raw, errors='coerce')
        an_val = pd.to_numeric(an_raw, errors='coerce')
        hom_val = pd.to_numeric(hom_raw, errors='coerce') if hom_raw is not None else None
        prepared.append((int(idx), variant_id, files_field, file_count, ac_val, an_val, hom_val))
    return prepared


# ----------------------------
# Main parallel processing + chunked flush
# ----------------------------
def run_per_variant_tests_parallel(df, config):
    print("Preparing data for parallel processing...")
    initialize_output_files(config)
    write_headers = True

    prepared_data = prepare_variant_data(df, config)
    total_variants = len(prepared_data)
    del df
    gc.collect()

    print(f"Processing {total_variants} variants using {config['n_workers']} threads...")
    print(f"Chunk size: {config['chunk_size']} variants per chunk")
    print(f"Barnard test settings:")
    print(f"  - Max total sample size for exact: {config.get('barnard_max_total', 10000)}")
    print(f"  - Downsample controls to: {config.get('barnard_downsample_to', 5000)}")
    print(f"  - Timeout: {config.get('barnard_timeout', 10)} seconds\n")

    process_func = partial(process_single_variant, config=config)

    total_fisher_count = 0
    total_barnard_count = 0
    total_bayes_count = 0
    total_skipped = 0

    chunk_size = int(config.get('chunk_size', 200))
    n_chunks = (total_variants + chunk_size - 1) // chunk_size

    with ThreadPoolExecutor(max_workers=config['n_workers']) as executor:
        with tqdm(total=n_chunks, desc="Processing chunks", unit="chunk") as pbar_chunks:
            for chunk_idx in range(n_chunks):
                start_idx = chunk_idx * chunk_size
                end_idx = min(start_idx + chunk_size, total_variants)
                chunk = prepared_data[start_idx:end_idx]

                fisher_rows = []
                barnard_rows = []
                bayes_rows = []
                skipped_rows = []

                with tqdm(total=len(chunk), desc=f" Chunk {chunk_idx+1}/{n_chunks}", unit="var", leave=False) as pbar_variants:
                    # map in parallel
                    for result in executor.map(process_func, chunk):
                        if result is None:
                            pbar_variants.update(1)
                            continue
                        if result.get('status') == 'processed':
                            fisher_rows.append(result['fisher'])
                            barnard_rows.append(result['barnard'])
                            bayes_rows.append(result['bayes'])
                        elif result.get('status') == 'skipped':
                            skipped_rows.append(result['skipped'])
                        pbar_variants.update(1)

                # convert to DataFrames
                df_fisher_chunk = pd.DataFrame(fisher_rows)
                df_barnard_chunk = pd.DataFrame(barnard_rows)
                df_bayes_chunk = pd.DataFrame(bayes_rows)
                df_skipped_chunk = pd.DataFrame(skipped_rows)

                total_fisher_count += len(df_fisher_chunk)
                total_barnard_count += len(df_barnard_chunk)
                total_bayes_count += len(df_bayes_chunk)
                total_skipped += len(df_skipped_chunk)

                # flush to disk
                if config.get('flush_every_chunk', True):
                    append_results_to_csv(df_fisher_chunk, config['output_variant_fisher'], write_headers, file_locks['fisher'])
                    append_results_to_csv(df_barnard_chunk, config['output_variant_barnard'], write_headers, file_locks['barnard'])
                    append_results_to_csv(df_bayes_chunk, config['output_variant_bayesian'], write_headers, file_locks['bayes'])
                    append_results_to_csv(df_skipped_chunk, config['output_skipped'], write_headers, file_locks['skipped'])
                    write_headers = False

                # free memory
                del fisher_rows, barnard_rows, bayes_rows, skipped_rows
                try:
                    del df_fisher_chunk, df_barnard_chunk, df_bayes_chunk, df_skipped_chunk
                except Exception:
                    pass
                gc.collect()

                pbar_chunks.update(1)

    print(f"\n✓ All chunks processed and flushed to disk.")
    print(f"Total written: Fisher={total_fisher_count}, Barnard={total_barnard_count}, Bayes={total_bayes_count}, Skipped={total_skipped}")

    # Apply multiple testing corrections by reading back files (if they exist)
    df_fisher_final = pd.DataFrame()
    df_barnard_final = pd.DataFrame()
    df_bayes_final = pd.DataFrame()

    # Fisher corrections
    if os.path.exists(config['output_variant_fisher']):
        df_fisher_final = pd.read_csv(config['output_variant_fisher'])
        if not df_fisher_final.empty:
            n_f = len(df_fisher_final)
            df_fisher_final['P_Value_Bonferroni'] = (df_fisher_final['P_Value'].fillna(1.0) * n_f).clip(upper=1.0)
            _, pvals_fdr, _, _ = multipletests(df_fisher_final['P_Value'].fillna(1.0), method='fdr_bh')
            df_fisher_final['P_Value_FDR'] = pvals_fdr
            df_fisher_final['Significant_FDR'] = df_fisher_final['P_Value_FDR'] < config['significance_level']
            df_fisher_final.to_csv(config['output_variant_fisher'], index=False)

    # Barnard corrections
    if os.path.exists(config['output_variant_barnard']):
        df_barnard_final = pd.read_csv(config['output_variant_barnard'])
        if not df_barnard_final.empty:
            n_b = len(df_barnard_final)
            df_barnard_final['P_Value_Bonferroni'] = (df_barnard_final['P_Value'].fillna(1.0) * n_b).clip(upper=1.0)
            _, pvals_fdr_b, _, _ = multipletests(df_barnard_final['P_Value'].fillna(1.0), method='fdr_bh')
            df_barnard_final['P_Value_FDR'] = pvals_fdr_b
            df_barnard_final['Significant_FDR'] = df_barnard_final['P_Value_FDR'] < config['significance_level']
            df_barnard_final.to_csv(config['output_variant_barnard'], index=False)
            
            # Report Barnard method statistics
            if 'Barnard_method' in df_barnard_final.columns:
                print("\nBarnard test methods used:")
                method_counts = df_barnard_final['Barnard_method'].value_counts()
                for method, count in method_counts.items():
                    print(f"  {method}: {count} variants")

    # Bayesian corrections
    if os.path.exists(config['output_variant_bayesian']):
        df_bayes_final = pd.read_csv(config['output_variant_bayesian'])
        if not df_bayes_final.empty:
            n_bv = len(df_bayes_final)
            df_bayes_final['P_Value_Bonferroni'] = (df_bayes_final['Bayes_p'].fillna(1.0) * n_bv).clip(upper=1.0)
            _, pvals_fdr_bv, _, _ = multipletests(df_bayes_final['Bayes_p'].fillna(1.0), method='fdr_bh')
            df_bayes_final['P_Value_FDR'] = pvals_fdr_bv
            df_bayes_final['Significant_FDR'] = df_bayes_final['P_Value_FDR'] < config['significance_level']
            df_bayes_final.to_csv(config['output_variant_bayesian'], index=False)

    gc.collect()
    return df_fisher_final, df_barnard_final, df_bayes_final


# ----------------------------
# Main entrypoint
# ----------------------------
def main():
    print("=" * 80)
    print("Per-variant strict pipeline (Homozygote-aware, multi-threaded)")
    print("Fixed Barnard implementation for large samples")
    print("=" * 80)
    print(f"Threads: {CONFIG['n_workers']}")
    print(f"Chunk size: {CONFIG['chunk_size']}")
    print(f"Flush per chunk: {CONFIG.get('flush_every_chunk', True)}")
    print(f"Input: {CONFIG['input_file']}")
    print()

    try:
        print("Loading input file...")
        df = pd.read_csv(CONFIG['input_file'], dtype=str, low_memory=False)
        print(f"Loaded {len(df)} variants")
    except Exception as e:
        print(f"Error loading input: {e}")
        return

    # Ensure required columns exist
    for col in ['SRR_Files', 'File_Count', 'Allele_Count',
                CONFIG['gnomad_ac_col'], CONFIG['gnomad_an_col'], CONFIG.get('gnomad_hom_col')]:
        if col not in df.columns:
            df[col] = np.nan

    df_fisher, df_barnard, df_bayes = run_per_variant_tests_parallel(df, CONFIG)

    print("\nFinal results:")
    if df_fisher is not None and not df_fisher.empty:
        print(f"✓ Fisher results rows: {len(df_fisher)} -> {CONFIG['output_variant_fisher']}")
    elif os.path.exists(CONFIG['output_variant_fisher']):
        print(f"✓ Fisher results written -> {CONFIG['output_variant_fisher']}")
    else:
        print("✗ No Fisher results (no variants with valid AC/AN)")

    if df_barnard is not None and not df_barnard.empty:
        print(f"✓ Barnard results rows: {len(df_barnard)} -> {CONFIG['output_variant_barnard']}")
    elif os.path.exists(CONFIG['output_variant_barnard']):
        print(f"✓ Barnard results written -> {CONFIG['output_variant_barnard']}")
    else:
        print("✗ No Barnard results (no variants with valid AC/AN or barnard missing)")

    if df_bayes is not None and not df_bayes.empty:
        print(f"✓ Bayesian results rows: {len(df_bayes)} -> {CONFIG['output_variant_bayesian']}")
    elif os.path.exists(CONFIG['output_variant_bayesian']):
        print(f"✓ Bayesian results written -> {CONFIG['output_variant_bayesian']}")
    else:
        print("✗ No Bayesian results (no variants with valid AC/AN)")

    # skipped file info
    if os.path.exists(CONFIG['output_skipped']):
        try:
            n_skip = len(pd.read_csv(CONFIG['output_skipped']))
            print(f"✓ Skipped variants log: {CONFIG['output_skipped']} (rows: {n_skip})")
        except Exception:
            print(f"✓ Skipped variants log exists: {CONFIG['output_skipped']}")
    else:
        print("✓ No skipped variants file created (no skips or file not written)")

    print("\n" + "=" * 80)


if __name__ == "__main__":
    main()
