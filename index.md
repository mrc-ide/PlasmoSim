# PlasmoSim

**PlasmoSim** is an R package for simulating epidemiological and genetic data from a lightweight model of *Plasmodium falciparum* transmission.

The package is designed for exploring how malaria transmission, mosquito ecology, human infection dynamics, parasite recombination, and movement between partially isolated populations shape patterns of parasite genetic structure.

## What does PlasmoSim do?

PlasmoSim simulates a simple but flexible malaria transmission system in which humans and mosquitoes interact across one or more demes: partially isolated sub-populations connected by migration. Within this framework, the package tracks both the epidemiological state of the system and the parasite genetic diversity generated through transmission and recombination.

The main simulation function, `sim_falciparum()`, models:

- human infection, recovery, and reinfection;
- mosquito infection, incubation, and survival;
- multiple simultaneous infections within hosts;
- transmission of parasite genotypes between mosquitoes and humans;
- meiotic recombination across multiple loci;
- migration between demes;
- scheduled sampling of individuals through time.

The model is intentionally lightweight. Rather than aiming to capture every biological detail of malaria transmission, PlasmoSim provides a fast and transparent simulation framework for asking spatial and population-genetic questions.

## Why use PlasmoSim?

PlasmoSim is useful when you want to explore questions such as:

- How does migration between populations affect parasite genetic relatedness?
- How quickly does genetic structure emerge between demes?
- How do mosquito population size, biting rate, and infection duration influence diversity?
- How does cotransmission of multiple parasite genotypes affect observed haplotypes?
- How do sampling time, sample size, and spatial sampling design affect genetic summaries?
- What patterns of relatedness might be expected under simple transmission scenarios?

Because the model is computationally efficient, it can be used to simulate multiple demes over long time periods, making it well suited to exploratory analyses, sensitivity analyses, and methodological testing.

## Simulation outputs

`sim_falciparum()` returns two main outputs.

The first is `daily_values`, a deme-level time series containing epidemiological summaries such as susceptible, exposed, and infectious humans and mosquitoes, together with the entomological inoculation rate.

The second is `indlevel`, an individual-level sampled dataset containing infection status and parasite haplotypes for sampled hosts. Haplotypes are represented across user-defined loci, allowing downstream analysis of genetic identity, relatedness, and structure.

## Genetic summaries

PlasmoSim includes helper functions for comparing simulated parasite haplotypes.

`get_haplotype_identity()` calculates the proportion of identical sites between two sets of haplotypes. This can be interpreted flexibly depending on how haplotypes are encoded. For example, it can represent average identity by state, or average identity by descent if loci encode ancestry.

`get_identity_matrix()` extends this idea to all sampled infections. It can return pairwise identity between individuals, or average identity within and between demes. This makes it straightforward to summarise spatial genetic structure from simulated data.

## A lightweight model for spatial malaria genetics

PlasmoSim is not intended to be a full malaria transmission model. Instead, it provides a deliberately simple simulation framework that links core epidemiological processes to parasite genetic outcomes.

This makes it particularly useful for developing intuition, testing analysis pipelines, and exploring how transmission assumptions might affect patterns seen in *P. falciparum* genomic data.
