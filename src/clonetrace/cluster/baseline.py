"""Variational Bayesian binomial mixture for single-sample VAF clustering."""

import math
from collections.abc import Sequence
from dataclasses import dataclass

from clonetrace.schemas import MutationRecord


@dataclass(frozen=True)
class MutationClusterAssignment:
    """Cluster assignment and posterior membership probabilities for one mutation."""

    mutation: MutationRecord
    cluster_id: int
    posterior_probability: float
    cluster_probabilities: tuple[float, ...] | None
    cellular_prevalence: float | None = None
    cellular_prevalence_std: float | None = None


@dataclass(frozen=True)
class BinomialCluster:
    """Posterior summary for one mixture component."""

    cluster_id: int
    posterior_mean_vaf: float
    mixture_weight: float
    mutation_count: int


@dataclass(frozen=True)
class BaselineClusteringResult:
    """Result from the Bayesian binomial mixture baseline."""

    sample: str
    method: str
    cluster_count: int
    assignments: tuple[MutationClusterAssignment, ...]
    clusters: tuple[BinomialCluster, ...]
    log_evidence_lower_bound: float
    iterations: int
    converged: bool


@dataclass(frozen=True)
class _VariationalFit:
    responsibilities: tuple[tuple[float, ...], ...]
    dirichlet_parameters: tuple[float, ...]
    beta_alpha_parameters: tuple[float, ...]
    beta_beta_parameters: tuple[float, ...]
    elbo: float
    iterations: int
    converged: bool


def _digamma(value: float) -> float:
    result = 0.0
    while value < 8:
        result -= 1 / value
        value += 1
    inverse = 1 / value
    inverse_squared = inverse * inverse
    result += (
        math.log(value)
        - inverse / 2
        - inverse_squared
        * (
            1 / 12
            - inverse_squared
            * (1 / 120 - inverse_squared * (1 / 252 - inverse_squared / 240))
        )
    )
    return result


def _log_beta(alpha: float, beta: float) -> float:
    return math.lgamma(alpha) + math.lgamma(beta) - math.lgamma(alpha + beta)


def _normalize_log_weights(log_weights: Sequence[float]) -> tuple[float, ...]:
    maximum = max(log_weights)
    weights = [math.exp(value - maximum) for value in log_weights]
    total = sum(weights)
    return tuple(weight / total for weight in weights)


def _initial_responsibilities(
    alt_counts: Sequence[int],
    ref_counts: Sequence[int],
    cluster_count: int,
) -> tuple[tuple[float, ...], ...]:
    observed_vafs = [alt / (alt + ref) for alt, ref in zip(alt_counts, ref_counts, strict=True)]
    sorted_vafs = sorted(observed_vafs)
    centers = [
        sorted_vafs[
            min(
                len(sorted_vafs) - 1,
                int((cluster + 0.5) * len(sorted_vafs) / cluster_count),
            )
        ]
        for cluster in range(cluster_count)
    ]
    epsilon = min(1e-3, 0.1 / cluster_count)
    rows: list[tuple[float, ...]] = []
    for vaf in observed_vafs:
        nearest = min(range(cluster_count), key=lambda cluster: abs(vaf - centers[cluster]))
        if cluster_count == 1:
            rows.append((1.0,))
        else:
            rows.append(
                tuple(
                    1 - epsilon * (cluster_count - 1) if cluster == nearest else epsilon
                    for cluster in range(cluster_count)
                )
            )
    return tuple(rows)


def _posterior_parameters(
    responsibilities: Sequence[Sequence[float]],
    alt_counts: Sequence[int],
    ref_counts: Sequence[int],
    *,
    cluster_count: int,
    dirichlet_concentration: float,
    beta_alpha: float,
    beta_beta: float,
) -> tuple[tuple[float, ...], tuple[float, ...], tuple[float, ...]]:
    effective_counts = [
        sum(row[cluster] for row in responsibilities) for cluster in range(cluster_count)
    ]
    alt_totals = [
        sum(row[cluster] * alt for row, alt in zip(responsibilities, alt_counts, strict=True))
        for cluster in range(cluster_count)
    ]
    ref_totals = [
        sum(row[cluster] * ref for row, ref in zip(responsibilities, ref_counts, strict=True))
        for cluster in range(cluster_count)
    ]
    return (
        tuple(count + dirichlet_concentration for count in effective_counts),
        tuple(count + beta_alpha for count in alt_totals),
        tuple(count + beta_beta for count in ref_totals),
    )


def _expected_logs(
    dirichlet_parameters: Sequence[float],
    beta_alpha_parameters: Sequence[float],
    beta_beta_parameters: Sequence[float],
) -> tuple[tuple[float, ...], tuple[float, ...], tuple[float, ...]]:
    total_dirichlet = sum(dirichlet_parameters)
    expected_log_weights = tuple(
        _digamma(value) - _digamma(total_dirichlet) for value in dirichlet_parameters
    )
    expected_log_vaf = tuple(
        _digamma(alpha) - _digamma(alpha + beta)
        for alpha, beta in zip(beta_alpha_parameters, beta_beta_parameters, strict=True)
    )
    expected_log_ref = tuple(
        _digamma(beta) - _digamma(alpha + beta)
        for alpha, beta in zip(beta_alpha_parameters, beta_beta_parameters, strict=True)
    )
    return expected_log_weights, expected_log_vaf, expected_log_ref


def _update_responsibilities(
    alt_counts: Sequence[int],
    ref_counts: Sequence[int],
    expected_log_weights: Sequence[float],
    expected_log_vaf: Sequence[float],
    expected_log_ref: Sequence[float],
) -> tuple[tuple[float, ...], ...]:
    rows: list[tuple[float, ...]] = []
    for alt, ref in zip(alt_counts, ref_counts, strict=True):
        log_probabilities = [
            expected_log_weights[cluster]
            + alt * expected_log_vaf[cluster]
            + ref * expected_log_ref[cluster]
            for cluster in range(len(expected_log_weights))
        ]
        rows.append(_normalize_log_weights(log_probabilities))
    return tuple(rows)


def _elbo(
    responsibilities: Sequence[Sequence[float]],
    alt_counts: Sequence[int],
    ref_counts: Sequence[int],
    dirichlet_parameters: Sequence[float],
    beta_alpha_parameters: Sequence[float],
    beta_beta_parameters: Sequence[float],
    *,
    dirichlet_concentration: float,
    beta_alpha: float,
    beta_beta: float,
) -> float:
    expected_log_weights, expected_log_vaf, expected_log_ref = _expected_logs(
        dirichlet_parameters,
        beta_alpha_parameters,
        beta_beta_parameters,
    )
    cluster_count = len(dirichlet_parameters)
    total_dirichlet = sum(dirichlet_parameters)
    value = 0.0
    for row_index, row in enumerate(responsibilities):
        alt = alt_counts[row_index]
        ref = ref_counts[row_index]
        depth = alt + ref
        log_binomial_coefficient = (
            math.lgamma(depth + 1) - math.lgamma(alt + 1) - math.lgamma(ref + 1)
        )
        for cluster, probability in enumerate(row):
            if probability == 0:
                continue
            expected_log_likelihood = (
                log_binomial_coefficient
                + alt * expected_log_vaf[cluster]
                + ref * expected_log_ref[cluster]
            )
            value += probability * (
                expected_log_likelihood + expected_log_weights[cluster] - math.log(probability)
            )

    value += (
        math.lgamma(cluster_count * dirichlet_concentration)
        - cluster_count * math.lgamma(dirichlet_concentration)
        - math.lgamma(total_dirichlet)
        + sum(math.lgamma(parameter) for parameter in dirichlet_parameters)
    )
    value += sum(
        (dirichlet_concentration - parameter) * expected_log_weight
        for parameter, expected_log_weight in zip(
            dirichlet_parameters,
            expected_log_weights,
            strict=True,
        )
    )

    for prior_alpha, prior_beta, posterior_alpha, posterior_beta, mean_log_vaf, mean_log_ref in zip(
        [beta_alpha] * cluster_count,
        [beta_beta] * cluster_count,
        beta_alpha_parameters,
        beta_beta_parameters,
        expected_log_vaf,
        expected_log_ref,
        strict=True,
    ):
        value += (
            _log_beta(posterior_alpha, posterior_beta)
            - _log_beta(prior_alpha, prior_beta)
            + (prior_alpha - posterior_alpha) * mean_log_vaf
            + (prior_beta - posterior_beta) * mean_log_ref
        )
    return value


def _fit_model(
    alt_counts: Sequence[int],
    ref_counts: Sequence[int],
    *,
    cluster_count: int,
    dirichlet_concentration: float,
    beta_alpha: float,
    beta_beta: float,
    max_iterations: int,
    tolerance: float,
) -> _VariationalFit:
    responsibilities = _initial_responsibilities(alt_counts, ref_counts, cluster_count)
    previous_elbo: float | None = None
    converged = False

    for _iteration in range(1, max_iterations + 1):
        dirichlet_parameters, beta_alpha_parameters, beta_beta_parameters = (
            _posterior_parameters(
                responsibilities,
                alt_counts,
                ref_counts,
                cluster_count=cluster_count,
                dirichlet_concentration=dirichlet_concentration,
                beta_alpha=beta_alpha,
                beta_beta=beta_beta,
            )
        )
        expected_logs = _expected_logs(
            dirichlet_parameters,
            beta_alpha_parameters,
            beta_beta_parameters,
        )
        updated_responsibilities = _update_responsibilities(
            alt_counts,
            ref_counts,
            *expected_logs,
        )
        dirichlet_parameters, beta_alpha_parameters, beta_beta_parameters = (
            _posterior_parameters(
                updated_responsibilities,
                alt_counts,
                ref_counts,
                cluster_count=cluster_count,
                dirichlet_concentration=dirichlet_concentration,
                beta_alpha=beta_alpha,
                beta_beta=beta_beta,
            )
        )
        current_elbo = _elbo(
            updated_responsibilities,
            alt_counts,
            ref_counts,
            dirichlet_parameters,
            beta_alpha_parameters,
            beta_beta_parameters,
            dirichlet_concentration=dirichlet_concentration,
            beta_alpha=beta_alpha,
            beta_beta=beta_beta,
        )
        responsibility_delta = max(
            abs(updated_responsibilities[row][cluster] - responsibilities[row][cluster])
            for row in range(len(responsibilities))
            for cluster in range(cluster_count)
        )
        elbo_delta = (
            abs(current_elbo - previous_elbo) if previous_elbo is not None else math.inf
        )
        responsibilities = updated_responsibilities
        if responsibility_delta <= tolerance and elbo_delta <= tolerance:
            converged = True
            break
        previous_elbo = current_elbo

    return _VariationalFit(
        responsibilities=responsibilities,
        dirichlet_parameters=dirichlet_parameters,
        beta_alpha_parameters=beta_alpha_parameters,
        beta_beta_parameters=beta_beta_parameters,
        elbo=current_elbo,
        iterations=_iteration,
        converged=converged,
    )


def fit_binomial_mixture(
    mutations: Sequence[MutationRecord],
    *,
    max_clusters: int = 6,
    dirichlet_concentration: float = 1.0,
    beta_alpha: float = 1.0,
    beta_beta: float = 1.0,
    max_iterations: int = 500,
    tolerance: float = 1e-6,
) -> BaselineClusteringResult:
    """Cluster one sample's mutations with a variational Bayesian binomial mixture.

    A symmetric Dirichlet prior regularizes mixture weights, and Beta priors
    regularize component VAFs. Candidate component counts from one through
    ``max_clusters`` are compared by their variational evidence lower bound.
    """

    if not mutations:
        raise ValueError("mutations must not be empty")
    if max_clusters < 1:
        raise ValueError("max_clusters must be at least 1")
    if max_iterations < 1:
        raise ValueError("max_iterations must be at least 1")
    if not math.isfinite(tolerance) or tolerance <= 0:
        raise ValueError("tolerance must be finite and greater than 0")
    for name, value in (
        ("dirichlet_concentration", dirichlet_concentration),
        ("beta_alpha", beta_alpha),
        ("beta_beta", beta_beta),
    ):
        if not math.isfinite(value) or value <= 0:
            raise ValueError(f"{name} must be finite and greater than 0")

    sample = mutations[0].sample
    if any(mutation.sample != sample for mutation in mutations):
        raise ValueError("all mutations must belong to the same sample")
    alt_counts = [mutation.alt_count for mutation in mutations]
    ref_counts = [mutation.ref_count for mutation in mutations]
    if any(alt + ref == 0 for alt, ref in zip(alt_counts, ref_counts, strict=True)):
        raise ValueError("mutations with zero total read depth cannot be clustered")

    candidate_count = min(max_clusters, len(mutations))
    fits = [
        _fit_model(
            alt_counts,
            ref_counts,
            cluster_count=cluster_count,
            dirichlet_concentration=dirichlet_concentration,
            beta_alpha=beta_alpha,
            beta_beta=beta_beta,
            max_iterations=max_iterations,
            tolerance=tolerance,
        )
        for cluster_count in range(1, candidate_count + 1)
    ]
    best_fit = max(fits, key=lambda fit: fit.elbo)
    old_cluster_count = len(best_fit.dirichlet_parameters)
    order = sorted(
        range(old_cluster_count),
        key=lambda cluster: (
            best_fit.beta_alpha_parameters[cluster]
            / (
                best_fit.beta_alpha_parameters[cluster]
                + best_fit.beta_beta_parameters[cluster]
            )
        ),
    )
    remapped_responsibilities = tuple(
        tuple(row[cluster] for cluster in order) for row in best_fit.responsibilities
    )
    posterior_mean_vafs = tuple(
        best_fit.beta_alpha_parameters[cluster]
        / (
            best_fit.beta_alpha_parameters[cluster]
            + best_fit.beta_beta_parameters[cluster]
        )
        for cluster in order
    )
    total_dirichlet = sum(best_fit.dirichlet_parameters)
    mixture_weights = tuple(
        best_fit.dirichlet_parameters[cluster] / total_dirichlet for cluster in order
    )
    assigned_counts = [0] * old_cluster_count
    assignments: list[MutationClusterAssignment] = []
    for mutation, probabilities in zip(mutations, remapped_responsibilities, strict=True):
        cluster_id = max(range(old_cluster_count), key=probabilities.__getitem__)
        assigned_counts[cluster_id] += 1
        assignments.append(
            MutationClusterAssignment(
                mutation=mutation,
                cluster_id=cluster_id,
                posterior_probability=probabilities[cluster_id],
                cluster_probabilities=probabilities,
            )
        )

    clusters = tuple(
        BinomialCluster(
            cluster_id=cluster_id,
            posterior_mean_vaf=posterior_mean_vafs[cluster_id],
            mixture_weight=mixture_weights[cluster_id],
            mutation_count=assigned_counts[cluster_id],
        )
        for cluster_id in range(old_cluster_count)
    )
    return BaselineClusteringResult(
        sample=sample,
        method="baseline",
        cluster_count=old_cluster_count,
        assignments=tuple(assignments),
        clusters=clusters,
        log_evidence_lower_bound=best_fit.elbo,
        iterations=best_fit.iterations,
        converged=best_fit.converged,
    )
