import pandas as pd


class ProcessRuninfoMetadata:
    @staticmethod
    def generate_runinfo_dict(runinfo_path, run_id_col='run_id'):
        run = pd.read_csv(runinfo_path)
        runinfo_rotated = pd.melt(run,
                                  id_vars=[run_id_col],
                                  value_vars=[colname for colname
                                              in list(run.columns)
                                              if colname != run_id_col],
                                  var_name='samples',
                                  value_name='value').dropna()

        runinfo_rotated['value'] = runinfo_rotated['value'].astype(int)

        runinfo_rotated['value'] = runinfo_rotated['value'].replace(
                                                                {1: 'case',
                                                                 0: 'control'})

        runinfo_rotated = runinfo_rotated.pivot(index=run_id_col,
                                                columns='value',
                                                values='samples')

        runinfo_rotated = runinfo_rotated.reset_index()

        runinfo_dict = runinfo_rotated.apply(lambda x:
                                             x.dropna().to_dict(),
                                             axis=1)

        return list(runinfo_dict)

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


if __name__ == '__main__':
    metadata = "/Users/pernillaericsson/Documents/medair1/apps/bio/dev_repos/parca/demo/runinfo/metadata.csv"
    runinfo = "/Users/pernillaericsson/Documents/medair1/apps/bio/dev_repos/parca/demo/runinfo/runinfo.csv"

    print(ProcessRuninfoMetadata.generate_runinfo_dict(runinfo))
    print(ProcessRuninfoMetadata.generate_metadata_dict(metadata))
