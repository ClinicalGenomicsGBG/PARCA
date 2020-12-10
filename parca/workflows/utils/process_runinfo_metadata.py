import pandas as pd


class ProcessRuninfoMetadata:
    @staticmethod
    def multiindex_pivot(df, index=None, columns=None, values=None):
        if index is None:
            names = list(df.index.names)
            df = df.reset_index()
        else:
            names = index
        list_index = df[names].values
        tuples_index = [tuple(i) for i in list_index]
        df = df.assign(tuples_index=tuples_index)
        df = df.pivot(index="tuples_index", columns=columns, values=values)
        tuples_index = df.index
        index = pd.MultiIndex.from_tuples(tuples_index, names=names)
        df.index = index
        return df

    @staticmethod
    def generate_runinfo_dict(runinfo_path, info_col=['run_id', 'start_date']):
        run = pd.read_csv(runinfo_path)

        runinfo_melt = pd.melt(run,
                               id_vars=info_col,
                               value_vars=[colname for colname
                                           in list(run.columns) if colname
                                           not in info_col],
                               var_name='samples',
                               value_name='value').dropna()
        runinfo_melt['value'] = runinfo_melt['value'].astype(int)
        runinfo_melt['value'] = runinfo_melt['value'].replace({1: 'case',
                                                               0: 'control'})

        runinfo_pivot = runinfo_melt.pipe(ProcessRuninfoMetadata.multiindex_pivot,
                                          index=info_col,
                                          columns='value', values='samples')
        runinfo_pivot = runinfo_pivot.reset_index()

        runinfo_dict = runinfo_pivot.apply(lambda x:
                                           x.dropna().to_dict(),
                                           axis=1)
        runinfo_dict = list(runinfo_dict)

        return runinfo_dict

    @staticmethod
    def generate_metadata_dict(metadata_path,
                               sample_id_col='sample_id',
                               fwd_or_rev_col='fwd_or_rev',
                               fwd_code="fwd",
                               rev_code="rev"):
        meta = pd.read_csv(metadata_path)

        detect_PE_or_SE = meta.groupby(sample_id_col).agg(
                        PE_or_SE=pd.NamedAgg(column=fwd_or_rev_col,
                                             aggfunc=lambda x: "PE"
                                             if(fwd_code in x.tolist() and
                                                rev_code in x.tolist())
                                             else "SE"))

        meta_PE_SE = pd.merge(meta, detect_PE_or_SE, on=[sample_id_col])

        meta_PE_SE = meta_PE_SE.fillna('NA')

        meta_PE_SE = meta_PE_SE.to_dict('records')

        return meta_PE_SE

    @staticmethod
    def nested_run_dict(list_run_dictionary,
                        start_date_col='start_date', run_id_col='run_id'):
        run_dict_nested = {}
        for run in list_run_dictionary:
            start_date_tmp = run[start_date_col]
            run_id_tmp = run[run_id_col]
            run_dict_nested[f'{start_date_tmp}_{run_id_tmp}'] = run
        return run_dict_nested

    @staticmethod
    def get_sample(run_dictionary, run_id, case_or_control):
        found_sample = run_dictionary.get(run_id).get(case_or_control)
        return found_sample

    @staticmethod
    def get_column(df, run_dictionary, run_id, case_or_control, column,
                   sample_id_col='sample_id', unique=False):
        found_column = df.loc[df.get(sample_id_col) == ProcessRuninfoMetadata.get_sample(run_dictionary, run_id, case_or_control)].get(column)
        if unique:
            found_column_list = found_column.unique().tolist()
            if len(found_column_list) != 1:
                return None
            else:
                return found_column_list[0]
        return found_column


if __name__ == '__main__':
    metadata = "/Users/pernillaericsson/Documents/medair1/apps/bio/dev_repos/parca/demo/runinfo/metadata.csv"
    runinfo = "/Users/pernillaericsson/Documents/medair1/apps/bio/dev_repos/parca/demo/runinfo/runinfo.csv"

    print(ProcessRuninfoMetadata.generate_runinfo_dict(runinfo))
    print(ProcessRuninfoMetadata.generate_metadata_dict(metadata))
