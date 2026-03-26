def run_single_sample(samples, params):
    results = []

    for sample, path in samples.items():
        adata = run_pipeline(sample, path, params)
        results.append(adata)

    return results