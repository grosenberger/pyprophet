import sys
import os
import click
import pandas as pd
import numpy as np
import sqlite3

from scipy.stats import rankdata
from .stats import error_statistics, lookup_values_from_error_table, final_err_table, summary_err_table
from .report import save_report
from shutil import copyfile
from .data_handling import check_sqlite_table


def statistics_report(data, outfile, context, analyte, parametric, pfdr, pi0_lambda, pi0_method, pi0_smooth_df, pi0_smooth_log_pi0, lfdr_truncate, lfdr_monotone, lfdr_transformation, lfdr_adj, lfdr_eps):

    error_stat, pi0 = error_statistics(data[data.decoy==0]['score'], data[data.decoy==1]['score'], parametric, pfdr, pi0_lambda, pi0_method, pi0_smooth_df, pi0_smooth_log_pi0, True, lfdr_truncate, lfdr_monotone, lfdr_transformation, lfdr_adj, lfdr_eps)

    stat_table = final_err_table(error_stat)
    summary_table = summary_err_table(error_stat)

    # print summary table
    click.echo("=" * 80)
    click.echo(summary_table)
    click.echo("=" * 80)

    p_values, s_values, peps, q_values = lookup_values_from_error_table(data["score"].values, error_stat)
    data["p_value"] = p_values
    data["s_value"] = s_values
    data["q_value"] = q_values
    data["pep"] = peps

    if context == 'run-specific':
        outfile = outfile + "_" + str(data['run_id'].unique()[0])

    # export PDF report
    save_report(outfile + "_" + context + "_" + analyte + ".pdf", outfile + ": " + context + " " + analyte + "-level error-rate control", data[data.decoy==1]["score"], data[data.decoy==0]["score"], stat_table["cutoff"], stat_table["svalue"], stat_table["qvalue"], data[data.decoy==0]["p_value"], pi0)

    return(data)


def infer_proteins(infile, outfile, context, parametric, pfdr, pi0_lambda, pi0_method, pi0_smooth_df, pi0_smooth_log_pi0, lfdr_truncate, lfdr_monotone, lfdr_transformation, lfdr_adj, lfdr_eps):

    con = sqlite3.connect(infile)

    if not check_sqlite_table(con, "SCORE_MS2"):
        raise click.ClickException("Apply scoring to MS2-level data before running protein-level scoring.")

    if context in ['global','experiment-wide','run-specific']:
        if context == 'global':
            run_id = 'NULL'
            group_id = 'PROTEIN.ID'
        else:
            run_id = 'RUN_ID'
            group_id = 'RUN_ID || "_" || PROTEIN.ID'

        con.executescript('''
CREATE INDEX IF NOT EXISTS idx_peptide_protein_mapping_protein_id ON PEPTIDE_PROTEIN_MAPPING (PROTEIN_ID);
CREATE INDEX IF NOT EXISTS idx_peptide_protein_mapping_peptide_id ON PEPTIDE_PROTEIN_MAPPING (PEPTIDE_ID);
CREATE INDEX IF NOT EXISTS idx_peptide_peptide_id ON PEPTIDE (ID);
CREATE INDEX IF NOT EXISTS idx_precursor_peptide_mapping_peptide_id ON PRECURSOR_PEPTIDE_MAPPING (PEPTIDE_ID);
CREATE INDEX IF NOT EXISTS idx_precursor_peptide_mapping_precursor_id ON PRECURSOR_PEPTIDE_MAPPING (PRECURSOR_ID);
CREATE INDEX IF NOT EXISTS idx_precursor_precursor_id ON PRECURSOR (ID);
CREATE INDEX IF NOT EXISTS idx_feature_precursor_id ON FEATURE (PRECURSOR_ID);
CREATE INDEX IF NOT EXISTS idx_feature_feature_id ON FEATURE (ID);
CREATE INDEX IF NOT EXISTS idx_score_ms2_feature_id ON SCORE_MS2 (FEATURE_ID);
''')

        data = pd.read_sql_query('''
SELECT %s AS RUN_ID,
       %s AS GROUP_ID,
       PROTEIN.ID AS PROTEIN_ID,
       PRECURSOR.DECOY AS DECOY,
       SCORE,
       "%s" AS CONTEXT
FROM PROTEIN
INNER JOIN
  (SELECT PEPTIDE_PROTEIN_MAPPING.PEPTIDE_ID AS PEPTIDE_ID,
          PROTEIN_ID
   FROM
     (SELECT PEPTIDE_ID,
             COUNT(*) AS NUM_PROTEINS
      FROM PEPTIDE_PROTEIN_MAPPING
      GROUP BY PEPTIDE_ID) AS PROTEINS_PER_PEPTIDE
   INNER JOIN PEPTIDE_PROTEIN_MAPPING ON PROTEINS_PER_PEPTIDE.PEPTIDE_ID = PEPTIDE_PROTEIN_MAPPING.PEPTIDE_ID
   WHERE NUM_PROTEINS == 1) AS PEPTIDE_PROTEIN_MAPPING ON PROTEIN.ID = PEPTIDE_PROTEIN_MAPPING.PROTEIN_ID
INNER JOIN PEPTIDE ON PEPTIDE_PROTEIN_MAPPING.PEPTIDE_ID = PEPTIDE.ID
INNER JOIN PRECURSOR_PEPTIDE_MAPPING ON PEPTIDE.ID = PRECURSOR_PEPTIDE_MAPPING.PEPTIDE_ID
INNER JOIN PRECURSOR ON PRECURSOR_PEPTIDE_MAPPING.PRECURSOR_ID = PRECURSOR.ID
INNER JOIN FEATURE ON PRECURSOR.ID = FEATURE.PRECURSOR_ID
INNER JOIN SCORE_MS2 ON FEATURE.ID = SCORE_MS2.FEATURE_ID
GROUP BY GROUP_ID
HAVING MAX(SCORE)
ORDER BY SCORE DESC
''' % (run_id, group_id, context), con)
    else:
        raise click.ClickException("Unspecified context selected.")

    data.columns = [col.lower() for col in data.columns]
    con.close()

    if context == 'run-specific':
        data = data.groupby('run_id').apply(statistics_report, outfile, context, "protein", parametric, pfdr, pi0_lambda, pi0_method, pi0_smooth_df, pi0_smooth_log_pi0, lfdr_truncate, lfdr_monotone, lfdr_transformation, lfdr_adj, lfdr_eps).reset_index()

    elif context in ['global', 'experiment-wide']:
        data = statistics_report(data, outfile, context, "protein", parametric, pfdr, pi0_lambda, pi0_method, pi0_smooth_df, pi0_smooth_log_pi0, lfdr_truncate, lfdr_monotone, lfdr_transformation, lfdr_adj, lfdr_eps)

    # store data in table
    if infile != outfile:
        copyfile(infile, outfile)

    con = sqlite3.connect(outfile)

    c = con.cursor()
    c.execute('SELECT count(name) FROM sqlite_master WHERE type="table" AND name="SCORE_PROTEIN"')
    if c.fetchone()[0] == 1:
        c.execute('DELETE FROM SCORE_PROTEIN WHERE CONTEXT =="%s"' % context)
    c.fetchall()

    df = data[['context','run_id','protein_id','score','p_value','q_value','pep']]
    df.columns = ['CONTEXT','RUN_ID','PROTEIN_ID','SCORE','PVALUE','QVALUE','PEP']
    table = "SCORE_PROTEIN"
    df.to_sql(table, con, index=False, dtype={"RUN_ID": "INTEGER"}, if_exists='append')

    con.close()


def infer_peptides(infile, outfile, context, parametric, pfdr, pi0_lambda, pi0_method, pi0_smooth_df, pi0_smooth_log_pi0, lfdr_truncate, lfdr_monotone, lfdr_transformation, lfdr_adj, lfdr_eps):

    con = sqlite3.connect(infile)

    if not check_sqlite_table(con, "SCORE_MS2"):
        raise click.ClickException("Apply scoring to MS2-level data before running peptide-level scoring.")

    if context in ['global','experiment-wide','run-specific']:
        if context == 'global':
            run_id = 'NULL'
            group_id = 'PEPTIDE.ID'
        else:
            run_id = 'RUN_ID'
            group_id = 'RUN_ID || "_" || PEPTIDE.ID'

        con.executescript('''
CREATE INDEX IF NOT EXISTS idx_peptide_peptide_id ON PEPTIDE (ID);
CREATE INDEX IF NOT EXISTS idx_precursor_peptide_mapping_peptide_id ON PRECURSOR_PEPTIDE_MAPPING (PEPTIDE_ID);
CREATE INDEX IF NOT EXISTS idx_precursor_peptide_mapping_precursor_id ON PRECURSOR_PEPTIDE_MAPPING (PRECURSOR_ID);
CREATE INDEX IF NOT EXISTS idx_precursor_precursor_id ON PRECURSOR (ID);
CREATE INDEX IF NOT EXISTS idx_feature_precursor_id ON FEATURE (PRECURSOR_ID);
CREATE INDEX IF NOT EXISTS idx_feature_feature_id ON FEATURE (ID);
CREATE INDEX IF NOT EXISTS idx_score_ms2_feature_id ON SCORE_MS2 (FEATURE_ID);
''')

        data = pd.read_sql_query('''
SELECT %s AS RUN_ID,
       %s AS GROUP_ID,
       PEPTIDE.ID AS PEPTIDE_ID,
       PRECURSOR.DECOY,
       SCORE,
       "%s" AS CONTEXT
FROM PEPTIDE
INNER JOIN PRECURSOR_PEPTIDE_MAPPING ON PEPTIDE.ID = PRECURSOR_PEPTIDE_MAPPING.PEPTIDE_ID
INNER JOIN PRECURSOR ON PRECURSOR_PEPTIDE_MAPPING.PRECURSOR_ID = PRECURSOR.ID
INNER JOIN FEATURE ON PRECURSOR.ID = FEATURE.PRECURSOR_ID
INNER JOIN SCORE_MS2 ON FEATURE.ID = SCORE_MS2.FEATURE_ID
GROUP BY GROUP_ID
HAVING MAX(SCORE)
ORDER BY SCORE DESC
''' % (run_id, group_id, context), con)
    else:
        raise click.ClickException("Unspecified context selected.")

    data.columns = [col.lower() for col in data.columns]
    con.close()

    if context == 'run-specific':
        data = data.groupby('run_id').apply(statistics_report, outfile, context, "peptide", parametric, pfdr, pi0_lambda, pi0_method, pi0_smooth_df, pi0_smooth_log_pi0, lfdr_truncate, lfdr_monotone, lfdr_transformation, lfdr_adj, lfdr_eps).reset_index()

    elif context in ['global', 'experiment-wide']:
        data = statistics_report(data, outfile, context, "peptide", parametric, pfdr, pi0_lambda, pi0_method, pi0_smooth_df, pi0_smooth_log_pi0, lfdr_truncate, lfdr_monotone, lfdr_transformation, lfdr_adj, lfdr_eps)

    # store data in table
    if infile != outfile:
        copyfile(infile, outfile)
    
    con = sqlite3.connect(outfile)

    c = con.cursor()
    c.execute('SELECT count(name) FROM sqlite_master WHERE type="table" AND name="SCORE_PEPTIDE"')
    if c.fetchone()[0] == 1:
        c.execute('DELETE FROM SCORE_PEPTIDE WHERE CONTEXT =="%s"' % context)
    c.fetchall()

    df = data[['context','run_id','peptide_id','score','p_value','q_value','pep']]
    df.columns = ['CONTEXT','RUN_ID','PEPTIDE_ID','SCORE','PVALUE','QVALUE','PEP']
    table = "SCORE_PEPTIDE"
    df.to_sql(table, con, index=False, dtype={"RUN_ID": "INTEGER"}, if_exists='append')

    con.close()


def infer_inter(infile, outfile, parametric, pfdr, pi0_lambda, pi0_method, pi0_smooth_df, pi0_smooth_log_pi0, lfdr_truncate, lfdr_monotone, lfdr_transformation, lfdr_adj, lfdr_eps):

    con = sqlite3.connect(infile)

    if not check_sqlite_table(con, "SCORE_MS2"):
        raise click.ClickException("Apply scoring to MS2-level data before running peptide-level scoring.")
    ipf_present = check_sqlite_table(con, "SCORE_IPF")

    con.executescript('''
CREATE INDEX IF NOT EXISTS idx_precursor_precursor_id ON PRECURSOR (ID);
CREATE INDEX IF NOT EXISTS idx_feature_precursor_id ON FEATURE (PRECURSOR_ID);
CREATE INDEX IF NOT EXISTS idx_feature_feature_id ON FEATURE (ID);
CREATE INDEX IF NOT EXISTS idx_score_ms2_feature_id ON SCORE_MS2 (FEATURE_ID);
CREATE INDEX IF NOT EXISTS idx_score_ipf_peptide_id ON SCORE_IPF (PEPTIDE_ID);
CREATE INDEX IF NOT EXISTS idx_score_ipf_feature_id ON SCORE_IPF (FEATURE_ID);
''')

    if ipf_present:
        click.echo("Info: Using IPF scores.")
        data = pd.read_sql_query('''
SELECT RUN_ID AS RUN_ID,
   SCORE_IPF.PEPTIDE_ID || '_' || PRECURSOR.ID AS PRECURSOR_ID,
   FEATURE.ID AS FEATURE_ID,
   FEATURE.DELTA_RT AS DELTA_RT,
   PRECURSOR.DECOY,
   SCORE_IPF.PEP
FROM SCORE_IPF
INNER JOIN FEATURE ON SCORE_IPF.FEATURE_ID = FEATURE.ID
INNER JOIN PRECURSOR ON FEATURE.PRECURSOR_ID = PRECURSOR.ID
INNER JOIN SCORE_MS2 ON FEATURE.ID = SCORE_MS2.FEATURE_ID
WHERE SCORE_MS2.RANK == 1 AND PRECURSOR.DECOY == 0
UNION
SELECT RUN_ID AS RUN_ID,
   PRECURSOR.ID AS PRECURSOR_ID,
   FEATURE.ID AS FEATURE_ID,
   FEATURE.DELTA_RT AS DELTA_RT,
   PRECURSOR.DECOY,
   PEP
FROM PRECURSOR
INNER JOIN FEATURE ON PRECURSOR.ID = FEATURE.PRECURSOR_ID
INNER JOIN SCORE_MS2 ON FEATURE.ID = SCORE_MS2.FEATURE_ID
WHERE SCORE_MS2.RANK == 1 AND PRECURSOR.DECOY == 1
''', con)

    else:
        data = pd.read_sql_query('''
SELECT RUN_ID AS RUN_ID,
   PRECURSOR.ID AS PRECURSOR_ID,
   FEATURE.ID AS FEATURE_ID,
   FEATURE.DELTA_RT AS DELTA_RT,
   PRECURSOR.DECOY,
   PEP
FROM PRECURSOR
INNER JOIN FEATURE ON PRECURSOR.ID = FEATURE.PRECURSOR_ID
INNER JOIN SCORE_MS2 ON FEATURE.ID = SCORE_MS2.FEATURE_ID
WHERE SCORE_MS2.RANK == 1
ORDER BY SCORE DESC
''', con)

    data.columns = [col.lower() for col in data.columns]
    con.close()

    # Select best scoring peptidoform per feature
    data = data.loc[data.groupby(['feature_id'])['pep'].idxmin()]

    # Compute inter-experiment alignment scores
    click.echo("Info: Compute Inter-Experiment Alignment (IEA) scores.")

    # Rank peak groups according expected retention time
    iea_ranked = data.copy()
    iea_ranked['delta_rt_abs'] = iea_ranked['delta_rt'].abs()
    iea_ranked['iearank'] = iea_ranked.groupby(['run_id','precursor_id'])['delta_rt'].rank(method='first').astype(int)

    iea_center = iea_ranked.loc[iea_ranked.groupby(['run_id','precursor_id'])['delta_rt_abs'].idxmin()][['run_id','precursor_id','iearank']]
    iea_center.columns = ['run_id','precursor_id','ieacenter']

    iea = iea_ranked.merge(iea_center, on=['run_id','precursor_id'])
    iea['rank_rt'] = iea['iearank'] - iea['ieacenter']

    # Compute IEA sum per precursor_id & rank_rt
    ieasum = iea.groupby(['precursor_id','rank_rt'])['pep'].apply(lambda x: sum(1-x)).reset_index(name='ieasum')
    iea = iea.merge(ieasum, on=['precursor_id','rank_rt'])

    # Compute IEA score per peak group
    iea['iea'] = iea.apply(lambda x: x['ieasum']-(1-x['pep']), axis=1)

    iea.loc[iea['iea'] > 15]['iea'] = 15
    iea.loc[iea['iea'] < -15]['iea'] = -15

    # Compute posterior error probabilities from IEA scores
    iea_error_stat, iea_pi0 = error_statistics(iea[iea.decoy==0]['iea'], iea[iea.decoy==1]['iea'], parametric, pfdr, pi0_lambda, pi0_method, pi0_smooth_df, pi0_smooth_log_pi0, True, lfdr_truncate, lfdr_monotone, lfdr_transformation, lfdr_adj, lfdr_eps)
    iea_p_values, _, iea_peps, _ = lookup_values_from_error_table(iea["iea"].values, iea_error_stat)
    iea_stat_table = final_err_table(iea_error_stat)

    iea["iea_p_value"] = iea_p_values
    iea["iea_pep"] = iea_peps

    # export IEA PDF report
    save_report(outfile + "_" + "inter_iea" + ".pdf", outfile + ": " + "IEA-level error-rate control", iea[iea.decoy==1]["iea"], iea[iea.decoy==0]["iea"], iea_stat_table["cutoff"], iea_stat_table["svalue"], iea_stat_table["qvalue"], iea[iea.decoy==0]["iea_p_value"], iea_pi0)

    # Reduce IEA scores
    iea = iea[['feature_id','iea_pep']]

    # Compute inter analyte & experiment probabilities
    click.echo("Info: Compute inter analyte & experiment probabilities.")
    inter = data.merge(iea, on='feature_id')

    # Bayesian integration
    # inter['inter'] = (((1-inter['iea_pep'])*(1-inter['nsi_pep'])*(1-inter['nsm_pep'])) * (1-inter['pep'])) / ((((1-inter['iea_pep'])*(1-inter['nsi_pep'])*(1-inter['nsm_pep'])) * (1-inter['pep'])) + ((inter['iea_pep']*inter['nsi_pep']*inter['nsm_pep']) * inter['pep']))
    inter['inter'] = ((1-inter['iea_pep']) * (1-inter['pep'])) / (((1-inter['iea_pep']) * (1-inter['pep'])) + (inter['iea_pep'] * inter['pep']))

    # Compute posterior error probabilities from inter scores
    inter_error_stat, inter_pi0 = error_statistics(inter[inter.decoy==0]['inter'], inter[inter.decoy==1]['inter'], parametric, pfdr, pi0_lambda, pi0_method, pi0_smooth_df, pi0_smooth_log_pi0, True, lfdr_truncate, lfdr_monotone, lfdr_transformation, lfdr_adj, lfdr_eps)
    inter_p_values, _, inter_peps, inter_q_values = lookup_values_from_error_table(inter["inter"].values, inter_error_stat)
    inter_stat_table = final_err_table(inter_error_stat)
    inter_summary_table = summary_err_table(inter_error_stat)

    inter["inter_p_value"] = inter_p_values
    inter["inter_pep"] = inter_peps
    inter["inter_q_value"] = inter_q_values

    # print summary table
    click.echo("=" * 80)
    click.echo(inter_summary_table)
    click.echo("=" * 80)

    # export inter PDF report
    save_report(outfile + "_" + "inter_integrated" + ".pdf", outfile + ": " + "inter-level error-rate control", inter[inter.decoy==1]["inter"], inter[inter.decoy==0]["inter"], inter_stat_table["cutoff"], inter_stat_table["svalue"], inter_stat_table["qvalue"], inter[inter.decoy==0]["inter_p_value"], inter_pi0)

    # store data in table
    if infile != outfile:
        copyfile(infile, outfile)
    
    con = sqlite3.connect(outfile)

    df = inter[['feature_id','inter_p_value','inter_q_value','inter_pep']]
    df.columns = ['FEATURE_ID','PVALUE','QVALUE','PEP']
    table = "SCORE_INTER"
    df.to_sql(table, con, index=False, dtype={"FEATURE_ID": "INTEGER"}, if_exists='replace')

    con.close()


def subsample_osw(infile, outfile, subsample_ratio, test):
    conn = sqlite3.connect(infile)
    ms1_present = check_sqlite_table(conn, "FEATURE_MS1")
    ms2_present = check_sqlite_table(conn, "FEATURE_MS2")
    transition_present = check_sqlite_table(conn, "FEATURE_TRANSITION")
    conn.close()

    conn = sqlite3.connect(outfile)
    c = conn.cursor()

    c.executescript('''
PRAGMA synchronous = OFF;

ATTACH DATABASE "%s" AS sdb;

CREATE TABLE RUN AS SELECT * FROM sdb.RUN;

DETACH DATABASE sdb;
''' % infile)
    click.echo("Info: Propagated runs of file %s to %s." % (infile, outfile))

    if subsample_ratio >= 1.0:
        c.executescript('''
ATTACH DATABASE "%s" AS sdb; 

CREATE TABLE FEATURE AS SELECT * FROM sdb.FEATURE; 

DETACH DATABASE sdb;
''' % infile)
    else:
        if test:
            c.executescript('''
ATTACH DATABASE "%s" AS sdb;

CREATE TABLE FEATURE AS 
SELECT *
FROM sdb.FEATURE
WHERE PRECURSOR_ID IN
    (SELECT ID
     FROM sdb.PRECURSOR
     LIMIT
       (SELECT ROUND(%s*COUNT(DISTINCT ID))
        FROM sdb.PRECURSOR));

DETACH DATABASE sdb;
''' % (infile, subsample_ratio))
        else:
            c.executescript('''
ATTACH DATABASE "%s" AS sdb;

CREATE TABLE FEATURE AS 
SELECT *
FROM sdb.FEATURE
WHERE PRECURSOR_ID IN
    (SELECT ID
     FROM sdb.PRECURSOR
     ORDER BY RANDOM()
     LIMIT
       (SELECT ROUND(%s*COUNT(DISTINCT ID))
        FROM sdb.PRECURSOR));

DETACH DATABASE sdb;
''' % (infile, subsample_ratio))
    click.echo("Info: Subsampled generic features of file %s to %s." % (infile, outfile))

    if ms1_present:
        if subsample_ratio >= 1.0:
            c.executescript('''
ATTACH DATABASE "%s" AS sdb;

CREATE TABLE FEATURE_MS1 AS 
SELECT *
FROM sdb.FEATURE_MS1;

DETACH DATABASE sdb;
''' % infile)
        else:
            c.executescript('''
ATTACH DATABASE "%s" AS sdb;

CREATE TABLE FEATURE_MS1 AS 
SELECT *
FROM sdb.FEATURE_MS1
WHERE sdb.FEATURE_MS1.FEATURE_ID IN
    (SELECT ID
     FROM FEATURE);

DETACH DATABASE sdb;
''' % infile)
        click.echo("Info: Subsampled MS1 features of file %s to %s." % (infile, outfile))

    if ms2_present:
        if subsample_ratio >= 1.0:
            c.executescript('''
ATTACH DATABASE "%s" AS sdb;

CREATE TABLE FEATURE_MS2 AS 
SELECT *
FROM sdb.FEATURE_MS2;

DETACH DATABASE sdb;
''' % infile)
        else:
            c.executescript('''
ATTACH DATABASE "%s" AS sdb;

CREATE TABLE FEATURE_MS2 AS 
SELECT *
FROM sdb.FEATURE_MS2
WHERE sdb.FEATURE_MS2.FEATURE_ID IN
    (SELECT ID
     FROM FEATURE);

DETACH DATABASE sdb;
''' % infile)
        click.echo("Info: Subsampled MS2 features of file %s to %s." % (infile, outfile))

    if transition_present:
        if subsample_ratio >= 1.0:
            c.executescript('''
ATTACH DATABASE "%s" AS sdb;

CREATE TABLE FEATURE_TRANSITION AS 
SELECT *
FROM sdb.FEATURE_TRANSITION;

DETACH DATABASE sdb;
''' % infile)
        else:
            c.executescript('''
ATTACH DATABASE "%s" AS sdb;

CREATE TABLE FEATURE_TRANSITION AS 
SELECT *
FROM sdb.FEATURE_TRANSITION
WHERE sdb.FEATURE_TRANSITION.FEATURE_ID IN
    (SELECT ID
     FROM FEATURE);

DETACH DATABASE sdb;
''' % infile)
        click.echo("Info: Subsampled transition features of file %s to %s." % (infile, outfile))

    conn.commit()
    conn.close()

    click.echo("Info: OSW file was subsampled.")


def reduce_osw(infile, outfile):
    conn = sqlite3.connect(infile)
    if not check_sqlite_table(conn, "SCORE_MS2"):
        raise click.ClickException("Apply scoring to MS2 data before reducing file for multi-run scoring.")
    ipf_present = check_sqlite_table(conn, "SCORE_IPF")
    conn.close()

    try:
        os.remove(outfile)
    except OSError:
        pass

    conn = sqlite3.connect(outfile)
    c = conn.cursor()

    query = '''
PRAGMA synchronous = OFF;

ATTACH DATABASE "%s" AS sdb;

CREATE TABLE RUN(ID INT PRIMARY KEY NOT NULL,
                 FILENAME TEXT NOT NULL);

INSERT INTO RUN
SELECT *
FROM sdb.RUN;

CREATE TABLE SCORE_MS2(FEATURE_ID INTEGER, SCORE REAL, RANK INTEGER, PEP REAL);

INSERT INTO SCORE_MS2 (FEATURE_ID, SCORE, RANK, PEP)
SELECT FEATURE_ID,
       SCORE,
       RANK,
       PEP
FROM sdb.SCORE_MS2
WHERE RANK == 1;

CREATE TABLE FEATURE(ID INT PRIMARY KEY NOT NULL,
                     RUN_ID INT NOT NULL,
                     PRECURSOR_ID INT NOT NULL,
                     DELTA_RT REAL);

INSERT INTO FEATURE (ID, RUN_ID, PRECURSOR_ID, DELTA_RT)
SELECT ID,
       RUN_ID,
       PRECURSOR_ID,
       DELTA_RT
FROM sdb.FEATURE
WHERE ID IN
    (SELECT FEATURE_ID
     FROM SCORE_MS2);
''' % infile

    if ipf_present:
        query += '''
CREATE TABLE SCORE_IPF(FEATURE_ID INTEGER, PEPTIDE_ID INTEGER, PEP REAL);

INSERT INTO SCORE_IPF (FEATURE_ID, PEPTIDE_ID, PEP)
SELECT FEATURE_ID,
       PEPTIDE_ID,
       PEP
FROM sdb.SCORE_IPF;
'''

    c.executescript(query)

    conn.commit()
    conn.close()

    click.echo("Info: OSW file was reduced for multi-run scoring.")


def merge_osw(infiles, outfile, templatefile, same_run):
    conn = sqlite3.connect(infiles[0])
    reduced = check_sqlite_table(conn, "SCORE_MS2")
    conn.close()

    if reduced:
        merge_oswr(infiles, outfile, templatefile, same_run)
    else:
        merge_osws(infiles, outfile, templatefile, same_run)


def merge_osws(infiles, outfile, templatefile, same_run):
    # Copy the first file to have a template
    copyfile(templatefile, outfile)
    conn = sqlite3.connect(outfile)
    c = conn.cursor()
    if same_run:
        c.execute("SELECT ID, FILENAME FROM RUN")
        result = c.fetchall()
        if len(result) != 1:
            raise click.ClickException("Input for same-run merge contains more than one run.")
        runid, rname = result[0]

    c.executescript('''
PRAGMA synchronous = OFF;

DROP TABLE IF EXISTS RUN;

DROP TABLE IF EXISTS FEATURE;

DROP TABLE IF EXISTS FEATURE_MS1;

DROP TABLE IF EXISTS FEATURE_MS2;

DROP TABLE IF EXISTS FEATURE_TRANSITION;

DROP TABLE IF EXISTS SCORE_MS1;

DROP TABLE IF EXISTS SCORE_MS2;

DROP TABLE IF EXISTS SCORE_TRANSITION;

DROP TABLE IF EXISTS SCORE_PEPTIDE;

DROP TABLE IF EXISTS SCORE_PROTEIN;

DROP TABLE IF EXISTS SCORE_IPF;

ATTACH DATABASE "%s" AS sdb;

CREATE TABLE RUN AS SELECT * FROM sdb.RUN LIMIT 0;

CREATE TABLE FEATURE AS SELECT * FROM sdb.FEATURE LIMIT 0;

CREATE TABLE FEATURE_MS1 AS SELECT * FROM sdb.FEATURE_MS1 LIMIT 0;

CREATE TABLE FEATURE_MS2 AS SELECT * FROM sdb.FEATURE_MS2 LIMIT 0;

CREATE TABLE FEATURE_TRANSITION AS SELECT * FROM sdb.FEATURE_TRANSITION LIMIT 0;

DETACH DATABASE sdb;
''' % (infiles[0]))

    for infile in infiles:
        # Only create a single run entry (all files are presumably from the same run)
        if same_run:
            c.executescript('''INSERT INTO RUN (ID, FILENAME) VALUES (%s, '%s')''' % (runid, rname) )
            break;
        else:
            c.executescript('''
ATTACH DATABASE "%s" AS sdb;

INSERT INTO RUN SELECT * FROM sdb.RUN;

DETACH DATABASE sdb;
''' % infile)

        click.echo("Info: Merged runs of file %s to %s." % (infile, outfile))

    # Now merge the run-specific data into the output file:
    #   Note: only tables FEATURE, FEATURE_MS1, FEATURE_MS2 and FEATURE_TRANSITION are run-specific
    for infile in infiles:
        c.executescript('''
ATTACH DATABASE "%s" AS sdb; 

INSERT INTO FEATURE SELECT * FROM sdb.FEATURE; 

DETACH DATABASE sdb;
''' % infile)
        click.echo("Info: Merged generic features of file %s to %s." % (infile, outfile))
    if same_run:
        # Fix run id assuming we only have a single run
        c.executescript('''UPDATE FEATURE SET RUN_ID = %s''' % runid)

    for infile in infiles:
        c.executescript('''
ATTACH DATABASE "%s" AS sdb;

INSERT INTO FEATURE_MS1
SELECT *
FROM sdb.FEATURE_MS1;

DETACH DATABASE sdb;
''' % infile)
        click.echo("Info: Merged MS1 features of file %s to %s." % (infile, outfile))

    for infile in infiles:
        c.executescript('''
ATTACH DATABASE "%s" AS sdb;

INSERT INTO FEATURE_MS2
SELECT *
FROM sdb.FEATURE_MS2;

DETACH DATABASE sdb;
''' % infile)
        click.echo("Info: Merged MS2 features of file %s to %s." % (infile, outfile))

    for infile in infiles:
        c.executescript('''
ATTACH DATABASE "%s" AS sdb;

INSERT INTO FEATURE_TRANSITION
SELECT *
FROM sdb.FEATURE_TRANSITION;

DETACH DATABASE sdb;
''' % infile)
        click.echo("Info: Merged transition features of file %s to %s." % (infile, outfile))

    conn.commit()
    conn.close()

    click.echo("Info: All OSWS files were merged.")


def merge_oswr(infiles, outfile, templatefile, same_run):
    # Copy the template to the output file
    copyfile(templatefile, outfile)
    conn = sqlite3.connect(outfile)
    c = conn.cursor()
    if same_run:
        c.execute("SELECT ID, FILENAME FROM RUN")
        result = c.fetchall()
        if len(result) != 1:
            raise click.ClickException("Input for same-run merge contains more than one run.")
        runid, rname = result[0]

    c.executescript('''
PRAGMA synchronous = OFF;

DROP TABLE IF EXISTS RUN;

DROP TABLE IF EXISTS FEATURE;

DROP TABLE IF EXISTS FEATURE_MS1;

DROP TABLE IF EXISTS FEATURE_MS2;

DROP TABLE IF EXISTS FEATURE_TRANSITION;

DROP TABLE IF EXISTS SCORE_MS1;

DROP TABLE IF EXISTS SCORE_MS2;

DROP TABLE IF EXISTS SCORE_TRANSITION;

DROP TABLE IF EXISTS SCORE_PEPTIDE;

DROP TABLE IF EXISTS SCORE_PROTEIN;

DROP TABLE IF EXISTS SCORE_IPF;

CREATE TABLE RUN(ID INT PRIMARY KEY NOT NULL,
                 FILENAME TEXT NOT NULL);

CREATE TABLE SCORE_MS2(FEATURE_ID INTEGER, SCORE REAL, RANK INTEGER, PEP REAL);

CREATE TABLE FEATURE(ID INT PRIMARY KEY NOT NULL,
                     RUN_ID INT NOT NULL,
                     PRECURSOR_ID INT NOT NULL,
                     DELTA_RT DOUBLE NOT NULL);

CREATE TABLE SCORE_IPF(FEATURE_ID INTEGER, PEPTIDE_ID INTEGER, PEP REAL);
''')

    for infile in infiles:
        # Only create a single run entry (all files are presumably from the same run)
        if same_run:
            c.executescript('''INSERT INTO RUN (ID, FILENAME) VALUES (%s, '%s')''' % (runid, rname) )
            break;
        else:
            c.executescript('ATTACH DATABASE "%s" AS sdb; INSERT INTO RUN SELECT * FROM sdb.RUN; DETACH DATABASE sdb;' % infile)
        click.echo("Info: Merged runs of file %s to %s." % (infile, outfile))

    for infile in infiles:
        c.executescript('ATTACH DATABASE "%s" AS sdb; INSERT INTO FEATURE SELECT * FROM sdb.FEATURE; DETACH DATABASE sdb;' % infile)
        click.echo("Info: Merged generic features of file %s to %s." % (infile, outfile))

    if same_run:
        # Fix run id assuming we only have a single run
        c.executescript('''UPDATE FEATURE SET RUN_ID = %s''' % runid)

    for infile in infiles:
        c.executescript('ATTACH DATABASE "%s" AS sdb; INSERT INTO SCORE_MS2 SELECT * FROM sdb.SCORE_MS2; DETACH DATABASE sdb;' % infile)
        click.echo("Info: Merged MS2 scores of file %s to %s." % (infile, outfile))

    for infile in infiles:
        conn2 = sqlite3.connect(infile)
        ipf_present = check_sqlite_table(conn2, "SCORE_IPF")
        conn2.close()

        if ipf_present:
            c.executescript('ATTACH DATABASE "%s" AS sdb; INSERT INTO SCORE_IPF SELECT * FROM sdb.SCORE_IPF; DETACH DATABASE sdb;' % infile)
            click.echo("Info: Merged IPF scores of file %s to %s." % (infile, outfile))

    # Drop IPF table if empty
    c.execute('SELECT count(*) FROM SCORE_IPF')
    if c.fetchone()[0] == 0:
        c.execute('DROP TABLE SCORE_IPF;')
    c.fetchall()

    conn.commit()
    conn.close()

    click.echo("Info: All reduced OSWR files were merged.")


def backpropagate_oswr(infile, outfile, apply_scores):
    # store data in table
    if infile != outfile:
        copyfile(infile, outfile)

    # find out what tables exist in the scores
    score_con = sqlite3.connect(apply_scores)
    inter_present = check_sqlite_table(score_con, "SCORE_INTER")
    peptide_present = check_sqlite_table(score_con, "SCORE_PEPTIDE")
    protein_present = check_sqlite_table(score_con, "SCORE_PROTEIN")
    score_con.close()
    if not (inter_present or peptide_present or protein_present):
        raise click.ClickException('Backpropagation requires inter, peptide or protein-level contexts.')

    # build up the list
    script = list()
    script.append('PRAGMA synchronous = OFF;')
    script.append('DROP TABLE IF EXISTS SCORE_INTER;')
    script.append('DROP TABLE IF EXISTS SCORE_PEPTIDE;')
    script.append('DROP TABLE IF EXISTS SCORE_PROTEIN;')

    # create the tables
    if inter_present:
        script.append('CREATE TABLE SCORE_INTER (FEATURE_ID INTEGER, PVALUE REAL, QVALUE REAL, PEP REAL);')
    if peptide_present:
        script.append('CREATE TABLE SCORE_PEPTIDE (CONTEXT TEXT, RUN_ID INTEGER, PEPTIDE_ID INTEGER, SCORE REAL, PVALUE REAL, QVALUE REAL, PEP REAL);')
    if protein_present:
        script.append('CREATE TABLE SCORE_PROTEIN (CONTEXT TEXT, RUN_ID INTEGER, PROTEIN_ID INTEGER, SCORE REAL, PVALUE REAL, QVALUE REAL, PEP REAL);')

    # copy across the tables
    script.append('ATTACH DATABASE "{}" AS sdb;'.format(apply_scores))
    insert_table_fmt = 'INSERT INTO {0}\nSELECT *\nFROM sdb.{0};'
    if inter_present:
        script.append(insert_table_fmt.format('SCORE_INTER'))
    if peptide_present:
        script.append(insert_table_fmt.format('SCORE_PEPTIDE'))
    if protein_present:
        script.append(insert_table_fmt.format('SCORE_PROTEIN'))

    # execute the script
    conn = sqlite3.connect(outfile)
    c = conn.cursor()
    c.executescript('\n'.join(script))
    conn.commit()
    conn.close()

    click.echo("Info: All multi-run data was backpropagated.")
