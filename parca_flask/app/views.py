#! /usr/bin/python3
# Maintainer Pernilla Ericsson

from flask import Flask, render_template, url_for, request, send_from_directory, abort

from app import app

import pandas as pd
import numpy as np
import glob
import os
import re


app.config["main_page_stats"] = "/data/main_page/main_page_stats_all.tsv" #"/Users/pernillaericsson/Documents/medair1/medstore/logs/pipeline_logfiles/parca/webinterface/main_page/main_page_stats_all.tsv"
app.config["run_id_dir"] = "/data" #"/Users/pernillaericsson/Documents/medair1/medstore/logs/pipeline_logfiles/parca/webinterface"

app.config["krona_expr"] = "_web/krona/*.krona.html"
app.config["tableview_expr"] = "_web/tableview/*_readcount_tableview.tsv"
app.config["detailed_stats_expr"] = "_web/tableview/*_detailed_stats.tsv"

app.config["dir_expr"] = "_web"
app.config["organism_fastq_expr"] = "organism_fastq"
app.config["kingdom_fastq_expr"] = "kingdom_fastq"
app.config["unclassified_fastq_expr"] = "unclassified_fastq/unclassified.fastq.gz"

app.config["fastq_suffix"] = ".fastq.gz"

@app.route('/')
def index():
    df = pd.read_csv(app.config["main_page_stats"], sep="\t")
    df = df.fillna("NA")

    return render_template("index.html",
                column_names=df.columns.values,
                row_lists=list(df.values.tolist()),
                zip=zip)

@app.route('/krona')
def krona():
    selected_sample = request.args.get('type')
    krona_path = glob.glob(os.path.join(app.config["run_id_dir"], f'{selected_sample}{app.config["krona_expr"]}'))[0]
    krona_relpath = re.sub(app.config["run_id_dir"], '', krona_path).strip('/')
    
    detailed_stats = glob.glob(os.path.join(app.config["run_id_dir"], f'{selected_sample}{app.config["detailed_stats_expr"]}'))[0]
    df = pd.read_csv(detailed_stats, sep="\t")

    options = ["count_raw_reads", "count_bbduk_trimmed_reads", "count_kmer_input", "count_unclassified_reads"]
    kmer_options_step = ["count_species_genus_higher"]
    kmer_options_type = ["species"]
    blast_options_step = ["count_reads_taxid_SubsetBLAST", "count_reads_taxid_BLASTnt"] 
    blast_options_type = ["reads"]

    df_subset1 = df[df['processing_step'].isin(options)]
    df_subset2 = df[df['processing_step'].isin(kmer_options_step) & df['type'].isin(kmer_options_type)]
    df_subset3 = df[df['processing_step'].isin(blast_options_step) & df['type'].isin(blast_options_type)]

    df_summary_full = pd.concat([df_subset1, df_subset2, df_subset3])
    df_summary_full['count'] = df_summary_full.drop(['processing_step', 'type'], axis=1).sum(axis=1).astype(int)
    df_summary = df_summary_full[['processing_step','count']]
    df_summary = df_summary.fillna("NA")
    df_summary.set_index('processing_step', inplace=True)
    selected_sample_record = df_summary.to_dict().get('count')

    return(render_template('krona.html', record = selected_sample_record, krona_path = krona_relpath, run_id = selected_sample) )

@app.route('/table_view')
def table_view():
    selected_sample = request.args.get('type')

    tableview_path = glob.glob(os.path.join(app.config["run_id_dir"], f'{selected_sample}{app.config["tableview_expr"]}'))[0]

    df = pd.read_csv(tableview_path, sep="\t")
    df = df.fillna("NA")

    all_column_names = df.columns
    
    df_group = df.groupby('superkingdom')
    df_group = sorted(df_group, key=lambda x: len(x[1]),reverse=True)
    df_lists = [list(sub_df.values.tolist()) for _,sub_df in df_group]


    detailed_stats = glob.glob(os.path.join(app.config["run_id_dir"], f'{selected_sample}{app.config["detailed_stats_expr"]}'))[0]
    df = pd.read_csv(detailed_stats, sep="\t")
    df = df.fillna("NA")

    case_sample = df[df['processing_step'] == 'case_sample'].get('type').item()
    case_path = os.path.join(app.config["run_id_dir"], f'{selected_sample}{app.config["dir_expr"]}', case_sample, app.config["organism_fastq_expr"])
    case_relpath = re.sub(app.config["run_id_dir"], '', case_path).strip('/')
    
    try:
        control_sample = df[df['processing_step'] == 'control_sample'].get('type').item()
        control_path = os.path.join(app.config["run_id_dir"], f'{selected_sample}{app.config["dir_expr"]}', control_sample, app.config["organism_fastq_expr"])
        control_relpath = re.sub(app.config["run_id_dir"], '', control_path).strip('/')
    except:
        control_relpath = None

    
    return(render_template('table_view.html', run_id = selected_sample,
           column_names=all_column_names,
           row_lists=df_lists,
           case=case_relpath,
           control=control_relpath,
           fastq_suffix=app.config["fastq_suffix"],
           zip=zip,
           len=len,
           str=str) ) 

@app.route('/detailed_stats')
def detailed_stats():
    # read_csv(f'detailed_stats_{selected_sample}')
    selected_sample = request.args.get('type')

    detailed_stats = glob.glob(os.path.join(app.config["run_id_dir"], f'{selected_sample}{app.config["detailed_stats_expr"]}'))[0]
    df = pd.read_csv(detailed_stats, sep="\t")
    df = df.fillna("NA")
    all_column_names = df.columns
    
    df_group = df.groupby('processing_step')
    df_lists = [list(sub_df.values.tolist()) for _,sub_df in df_group]

    case_sample = df[df['processing_step'] == 'case_sample'].get('type').item()

    case_path = os.path.join(app.config["run_id_dir"], f'{selected_sample}{app.config["dir_expr"]}', case_sample, app.config["kingdom_fastq_expr"])
    case_relpath = re.sub(app.config["run_id_dir"], '', case_path).strip('/')

    case_unclassified_path = os.path.join(app.config["run_id_dir"], f'{selected_sample}{app.config["dir_expr"]}', case_sample, app.config["unclassified_fastq_expr"])
    case_unclassified_relpath = re.sub(app.config["run_id_dir"], '', case_unclassified_path).strip('/')
    
    try:
        control_sample = df[df['processing_step'] == 'control_sample'].get('type').item()

        control_path = os.path.join(app.config["run_id_dir"], f'{selected_sample}{app.config["dir_expr"]}', control_sample, app.config["kingdom_fastq_expr"])
        control_relpath = re.sub(app.config["run_id_dir"], '', control_path).strip('/')

        control_unclassified_path = os.path.join(app.config["run_id_dir"], f'{selected_sample}{app.config["dir_expr"]}', control_sample, app.config["unclassified_fastq_expr"])
        control_unclassified_relpath = re.sub(app.config["run_id_dir"], '', control_unclassified_path).strip('/')
    except:
        control_relpath = None
        control_unclassified_relpath = None

    return(render_template('detailed_stats.html', run_id = selected_sample, 
           column_names=all_column_names,
           row_lists=df_lists,
           case=case_relpath,
           case_unclassified=case_unclassified_relpath,
           control=control_relpath,
           control_unclassified=control_unclassified_relpath,
           fastq_suffix=app.config["fastq_suffix"],
           len=len,
           zip=zip,
           int=int,
           str=str ) )

@app.route('/get_file/<path:file_name>')
def get_file(file_name):
    try:
        return send_from_directory(app.config["run_id_dir"], filename=file_name, as_attachment=False)
    except FileNotFoundError:
        abort(404)

@app.route('/get_file/<path:file_name>')
def download_file(file_name):
    try:
        return send_from_directory(app.config["run_id_dir"], filename=file_name, as_attachment=True)
    except FileNotFoundError:
        abort(404)

# if __name__ == '__main__':
#     app.run(debug=True, host='0.0.0.0')
