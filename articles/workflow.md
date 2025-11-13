# TTE simulation data manipulations

## Overview

We attempt to provide a big picture view of what is involved in a
clinical trial simulation with individual time-to-event data generated
for each patient. Primary interest is in group sequential trials,
usually with a single endpoint. However extensions could be made.

## Results data table

At the time of simulation planning an analysis plan for the trial is
needed. For a group sequential design, a data table to store results is
generated. Generally, the dimensions and variables planned for storage
would be planned up front. As a simple example, if there is a group
sequential design with 3 analyses planned, 15 data items for each
analysis and 10,000 simulations planned, a data table with 30,000 rows
and 15 columns has been used to store summary results. As each trial
simulation proceeds, a row is updated with results for each analysis.

## Simulated trial dataset generation

For each simulated trial, an initial table is generated with information
at a patient level. If trials are generated sequentially, the space
needed for this data table could be re-used, never requiring allocation
of more space. Each row contains data for a single patient. As an
example, we could simulate a trial with 500 patients and 10 data items
per patients. The data items would be in columns, the patients in rows.

## Dataset manipulations for analysis

Simulated trial data need to be manipulated to do any individual
analysis (interim or final) for a clinical trial. The following
operations are needed:

1.  Ordering data
2.  Selecting a subset for analysis
3.  Calculating individual patient results for the subset at the time of
    analysis.
4.  Performing statistical tests and computing treatment effect
    estimates as well as other summaries that will included in the
    results summary dataset described above. Types of computations
    needed are
    1.  Number of subjects by treatment group
    2.  Number of events by treatment group
    3.  Kaplan-Meier estimation of survival curves
    4.  Observed minus expected computations as well as weighting for
        logrank, weighted logrank calculations.
    5.  Using the survival package to compute hazard ratio estimates.

## Flow for simulating group sequential: one scenario algorithm

Group sequential design simulation flow:

- Generate a trial.
- Analyze repeatedly.
- Summarize across simulated trials.

![](figures/workflow.png)
