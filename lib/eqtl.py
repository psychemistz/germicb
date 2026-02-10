"""
eQTL reference data loaders: DICE, OneK1K, cohort eQTLs.
"""

from pathlib import Path
from typing import Optional

import numpy as np
import pandas as pd

from .config import DICE_HG38_DIR, ONEK1K_PATH, WES_DIR, FDR_THRESHOLD
from .log import log


def load_dice(dice_dir: Optional[Path] = None) -> Optional[pd.DataFrame]:
    """
    Load DICE immune cell eQTLs (hg38 lifted coordinates).

    Returns DataFrame with columns:
        chrom, pos, rsid, ref, alt, gene_symbol, gene_id, beta, pval, dice_celltype
    """
    dice_dir = dice_dir or DICE_HG38_DIR

    if not dice_dir.exists():
        log(f"  DICE hg38 directory not found: {dice_dir}")
        return None

    log(f"Loading DICE eQTLs from: {dice_dir}")

    files_to_load = list(dice_dir.glob('*_hg38.tsv'))
    if not files_to_load:
        log("  No DICE hg38 TSV files found")
        return None

    log(f"  Found {len(files_to_load)} cell type files")

    all_records = []

    for file_path in files_to_load:
        if file_path.stat().st_size == 0:
            continue

        dice_celltype = file_path.stem.replace('_hg38', '')
        df = pd.read_csv(file_path, sep='\t')

        for _, row in df.iterrows():
            info_dict = {}
            if pd.notna(row.get('INFO', '')):
                for item in str(row['INFO']).split(';'):
                    if '=' in item:
                        key, val = item.split('=', 1)
                        info_dict[key] = val

            if pd.notna(row.get('CHROM_hg38')) and pd.notna(row.get('POS_hg38')):
                all_records.append({
                    'chrom': row['CHROM_hg38'],
                    'pos': int(row['POS_hg38']),
                    'rsid': row.get('ID', ''),
                    'ref': row.get('REF', ''),
                    'alt': row.get('ALT', ''),
                    'gene_symbol': info_dict.get('GeneSymbol', ''),
                    'gene_id': info_dict.get('Gene', ''),
                    'beta': float(info_dict.get('Beta', 0)),
                    'pval': float(info_dict.get('Pvalue', 1)),
                    'dice_celltype': dice_celltype,
                })

    if not all_records:
        log("  No DICE records loaded")
        return None

    dice_df = pd.DataFrame(all_records)
    log(f"  Total DICE eQTLs: {len(dice_df):,} from {dice_df['dice_celltype'].nunique()} cell types")

    return dice_df


def load_onek1k(path: Optional[Path] = None, fdr_threshold: float = FDR_THRESHOLD) -> Optional[pd.DataFrame]:
    """
    Load OneK1K single-cell PBMC eQTLs (hg38).

    Source: Yazar et al., Science 2022 — ~1.27M PBMCs from 982 donors, 14 cell types.
    Independent of AlphaGenome training data.

    Returns DataFrame with columns:
        chrom, pos, rsid, gene, beta, pval, onek1k_celltype
    """
    path = path or ONEK1K_PATH

    if not path.exists():
        log(f"  OneK1K file not found: {path}")
        return None

    log(f"Loading OneK1K eQTLs: {path}")

    onek1k_df = pd.read_csv(path, sep='\t')
    log(f"  Loaded {len(onek1k_df):,} records")

    if 'FDR' in onek1k_df.columns:
        onek1k_df = onek1k_df[onek1k_df['FDR'] < fdr_threshold].copy()
        log(f"  After FDR < {fdr_threshold} filter: {len(onek1k_df):,}")

    onek1k_df['chrom'] = 'chr' + onek1k_df['CHR'].astype(str)
    onek1k_df['pos'] = onek1k_df['POS'].astype(int)
    onek1k_df['rsid'] = onek1k_df['RSID']
    onek1k_df['gene'] = onek1k_df['GENE']
    onek1k_df['beta'] = onek1k_df['SPEARMANS_RHO'].astype(float)
    onek1k_df['pval'] = onek1k_df['P_VALUE'].astype(float)
    onek1k_df['onek1k_celltype'] = onek1k_df['CELL_TYPE']

    log(f"  Cell types: {onek1k_df['onek1k_celltype'].nunique()}")
    log(f"  Unique variants: {onek1k_df['rsid'].nunique():,}")

    return onek1k_df[['chrom', 'pos', 'rsid', 'gene', 'beta', 'pval', 'onek1k_celltype']].copy()


def load_cohort_eqtls(
    wes_dir: Optional[Path] = None,
    fdr_threshold: float = FDR_THRESHOLD,
) -> Optional[pd.DataFrame]:
    """
    Load per-cohort local eQTL results for sanity check.

    These are from the WES immunotherapy cohort analyses (hg19 coordinates).
    Requires liftover to hg38 for matching against DICE/OneK1K/AlphaGenome.

    Returns DataFrame with columns:
        chrom, pos_hg19, pos_hg38, rsid, gene, beta, pval, fdr, cohort
    """
    wes_dir = wes_dir or WES_DIR

    if not wes_dir.exists():
        log(f"  WES directory not found: {wes_dir}")
        return None

    log(f"Loading cohort eQTLs from: {wes_dir}")

    try:
        from pyliftover import LiftOver
        lo = LiftOver('hg19', 'hg38')
    except ImportError:
        log("  WARNING: pyliftover not available, using rsid matching only")
        lo = None

    all_cohort_eqtls = []
    cohort_dirs = [d for d in wes_dir.iterdir()
                   if d.is_dir() and not d.name.startswith(('summary', 'logis', 'ridge', 'Trans'))]

    for cohort_dir in sorted(cohort_dirs):
        eqtl_file = cohort_dir / '9eQTL' / 'output' / 'fin_eqtl_tbl.txt'
        snp_annot_file = cohort_dir / '9eQTL' / 'data' / 'snp_annot.txt'

        if not eqtl_file.exists():
            continue

        try:
            eqtl_df = pd.read_csv(eqtl_file, sep=r'\s+')
        except Exception as e:
            log(f"  WARNING: Could not load {eqtl_file}: {e}")
            continue

        if eqtl_df.empty:
            continue

        eqtl_df = eqtl_df[eqtl_df['FDR'] < fdr_threshold].copy()
        if eqtl_df.empty:
            continue

        # Clean beta: replace Inf/-Inf with NaN
        eqtl_df['beta'] = pd.to_numeric(eqtl_df['beta'], errors='coerce')
        eqtl_df = eqtl_df[eqtl_df['beta'].notna() & np.isfinite(eqtl_df['beta'])].copy()

        eqtl_df['cohort'] = cohort_dir.name

        if snp_annot_file.exists():
            try:
                snp_annot = pd.read_csv(snp_annot_file, sep='\t')
                snp_annot = snp_annot.rename(columns={'CHROM': 'chrom', 'POS': 'pos_hg19', 'ID': 'rsid_annot'})
                eqtl_df = eqtl_df.merge(
                    snp_annot[['chrom', 'pos_hg19', 'rsid_annot']].drop_duplicates('rsid_annot'),
                    left_on='snps', right_on='rsid_annot', how='inner'
                )
            except Exception as e:
                log(f"  WARNING: Could not join snp_annot for {cohort_dir.name}: {e}")
                continue
        else:
            continue

        if eqtl_df.empty:
            continue

        # Liftover hg19 -> hg38
        if lo is not None:
            hg38_positions = []
            for _, row in eqtl_df.iterrows():
                chrom = row['chrom']
                pos = int(row['pos_hg19'])
                try:
                    result = lo.convert_coordinate(chrom, pos)
                    if result and len(result) > 0:
                        hg38_positions.append(int(result[0][1]))
                    else:
                        hg38_positions.append(np.nan)
                except Exception:
                    hg38_positions.append(np.nan)
            eqtl_df['pos_hg38'] = hg38_positions
            eqtl_df = eqtl_df[eqtl_df['pos_hg38'].notna()].copy()
            eqtl_df['pos_hg38'] = eqtl_df['pos_hg38'].astype(int)

        n_eqtls = len(eqtl_df)
        if n_eqtls > 0:
            log(f"  {cohort_dir.name}: {n_eqtls:,} significant eQTLs (FDR<{fdr_threshold})")
            all_cohort_eqtls.append(eqtl_df)

    if not all_cohort_eqtls:
        log("  No cohort eQTLs loaded")
        return None

    combined = pd.concat(all_cohort_eqtls, ignore_index=True)
    log(f"  Total cohort eQTLs: {len(combined):,} from {combined['cohort'].nunique()} cohorts")

    return combined
