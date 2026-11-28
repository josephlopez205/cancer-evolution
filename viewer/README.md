# CloneTrace Simulation Explorer

A small React + TypeScript viewer for the Phase 2 simulator output. The bundled sample is
simulated ground truth and is not an inferred result or a clinical/diagnostic output.

## Run locally

From this directory:

```powershell
npm install
npm run dev
```

The viewer loads `public/sample-bundle/ground_truth.json`,
`public/sample-bundle/mutations.tsv`, and
`public/sample-bundle/copy_number_segments.tsv`. The bundle version and core data
consistency are checked when the app loads; missing or invalid data shows an error state.

## Validate

```powershell
npm test
npm run lint
npm run build
```

This preview currently shows the simulator's tree, mutation read counts, and copy-number
events. Cohort, inferred cluster, benchmark-result, and timing views will be added only when
their pipeline outputs and versioned export contract are available.
