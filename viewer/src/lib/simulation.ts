export type Topology = 'linear' | 'branching'

export interface SimulationConfig {
  seed: number
  depth: number
  purity: number
  clone_count: number
  topology: Topology
  mutation_count: number
  cn_event_probability: number
  sample: string
}

export interface MutationObservation {
  sample: string
  gene: string
  chrom: string
  pos: number
  ref: string
  alt: string
  ref_count: number
  alt_count: number
}

export interface CopyNumberEvent {
  sample: string
  chrom: string
  start: number
  end: number
  total_cn: number
  minor_cn: number
}

export interface MutationAssignment {
  mutation_id: number
  clone: number
  gene: string
  chrom: string
  pos: number
  true_ccf: number
}

export interface SimulationBundle {
  config: SimulationConfig
  parents: number[]
  cloneFractions: number[]
  assignments: MutationAssignment[]
  mutations: MutationObservation[]
  segments: CopyNumberEvent[]
}

function isRecord(value: unknown): value is Record<string, unknown> {
  return typeof value === 'object' && value !== null && !Array.isArray(value)
}

function finiteNumber(value: unknown, label: string): number {
  if (typeof value !== 'number' || !Number.isFinite(value)) {
    throw new Error(`Invalid simulation data: ${label} must be a finite number.`)
  }
  return value
}

function parseTsv(text: string, label: string): Record<string, string>[] {
  const [headerLine, ...lines] = text.trimEnd().split(/\r?\n/)
  if (!headerLine) {
    throw new Error(`Invalid simulation data: ${label} has no header.`)
  }
  const headers = headerLine.split('\t')
  if (headers.some((header) => !header)) {
    throw new Error(`Invalid simulation data: ${label} contains an empty column name.`)
  }
  return lines.filter(Boolean).map((line, index) => {
    const values = line.split('\t')
    if (values.length !== headers.length) {
      throw new Error(`Invalid simulation data: malformed ${label} row ${index + 2}.`)
    }
    return Object.fromEntries(headers.map((header, valueIndex) => [header, values[valueIndex]]))
  })
}

function requiredText(row: Record<string, string>, field: string, label: string): string {
  const value = row[field]
  if (!value) {
    throw new Error(`Invalid simulation data: ${label}.${field} is missing.`)
  }
  return value
}

function parseInteger(row: Record<string, string>, field: string, label: string): number {
  const value = Number(requiredText(row, field, label))
  if (!Number.isSafeInteger(value) || value < 0) {
    throw new Error(`Invalid simulation data: ${label}.${field} must be a non-negative integer.`)
  }
  return value
}

export function parseSimulationBundle(
  groundTruthValue: unknown,
  mutationText: string,
  segmentText: string,
): SimulationBundle {
  if (!isRecord(groundTruthValue) || groundTruthValue.format_version !== '1') {
    throw new Error('Unsupported or invalid simulation bundle version.')
  }

  const configValue = groundTruthValue.config
  const treeValue = groundTruthValue.clone_tree
  const assignmentsValue = groundTruthValue.mutation_assignments
  if (!isRecord(configValue) || !isRecord(treeValue) || !Array.isArray(assignmentsValue)) {
    throw new Error('Invalid simulation data: required ground-truth fields are missing.')
  }

  const topology = configValue.topology
  if (topology !== 'linear' && topology !== 'branching') {
    throw new Error('Invalid simulation data: unsupported topology.')
  }

  const config: SimulationConfig = {
    seed: finiteNumber(configValue.seed, 'config.seed'),
    depth: finiteNumber(configValue.depth, 'config.depth'),
    purity: finiteNumber(configValue.purity, 'config.purity'),
    clone_count: finiteNumber(configValue.clone_count, 'config.clone_count'),
    topology,
    mutation_count: finiteNumber(configValue.mutation_count, 'config.mutation_count'),
    cn_event_probability: finiteNumber(
      configValue.cn_event_probability,
      'config.cn_event_probability',
    ),
    sample: typeof configValue.sample === 'string' ? configValue.sample : '',
  }
  if (!config.sample || config.clone_count < 1 || config.depth < 1 || config.mutation_count < 1) {
    throw new Error('Invalid simulation data: configuration values are out of range.')
  }

  const parents = treeValue.parents
  const fractions = treeValue.clone_fractions
  if (
    !Array.isArray(parents) ||
    !Array.isArray(fractions) ||
    parents.length !== config.clone_count ||
    fractions.length !== config.clone_count
  ) {
    throw new Error('Invalid simulation data: clone tree does not match the configuration.')
  }

  const normalizedParents = parents.map((parent, index) => {
    if (!Number.isSafeInteger(parent) || parent < 0 || parent >= parents.length) {
      throw new Error(`Invalid simulation data: clone ${index} has an invalid parent.`)
    }
    if (index === 0 ? parent !== 0 : parent >= index) {
      throw new Error(`Invalid simulation data: clone ${index} is not rooted correctly.`)
    }
    return parent
  })
  const cloneFractions = fractions.map((fraction, index) => {
    const value = finiteNumber(fraction, `clone_tree.clone_fractions[${index}]`)
    if (value < 0 || value > 1) {
      throw new Error('Invalid simulation data: clone fractions must be between 0 and 1.')
    }
    return value
  })
  if (Math.abs(cloneFractions.reduce((sum, value) => sum + value, 0) - 1) > 1e-6) {
    throw new Error('Invalid simulation data: clone fractions must sum to 1.')
  }

  const mutations = parseTsv(mutationText, 'mutations.tsv').map((row, index) => ({
    sample: requiredText(row, 'sample', `mutations.tsv row ${index + 2}`),
    gene: requiredText(row, 'gene', `mutations.tsv row ${index + 2}`),
    chrom: requiredText(row, 'chrom', `mutations.tsv row ${index + 2}`),
    pos: parseInteger(row, 'pos', `mutations.tsv row ${index + 2}`),
    ref: requiredText(row, 'ref', `mutations.tsv row ${index + 2}`),
    alt: requiredText(row, 'alt', `mutations.tsv row ${index + 2}`),
    ref_count: parseInteger(row, 'ref_count', `mutations.tsv row ${index + 2}`),
    alt_count: parseInteger(row, 'alt_count', `mutations.tsv row ${index + 2}`),
  }))
  if (mutations.some((mutation) => mutation.sample !== config.sample || mutation.pos < 1)) {
    throw new Error('Invalid simulation data: mutation rows do not match the sample configuration.')
  }
  const segments = parseTsv(segmentText, 'copy_number_segments.tsv').map((row, index) => {
    const label = `copy_number_segments.tsv row ${index + 2}`
    const totalCn = Number(requiredText(row, 'total_cn', label))
    const minorCn = Number(requiredText(row, 'minor_cn', label))
    const start = parseInteger(row, 'start', label)
    const end = parseInteger(row, 'end', label)
    if (
      !Number.isFinite(totalCn) ||
      !Number.isFinite(minorCn) ||
      totalCn < 0 ||
      minorCn < 0 ||
      minorCn > totalCn ||
      end <= start
    ) {
      throw new Error(`Invalid simulation data: ${label} has invalid copy-number values.`)
    }
    if (requiredText(row, 'sample', label) !== config.sample) {
      throw new Error(`Invalid simulation data: ${label} does not match the sample configuration.`)
    }
    return {
      sample: config.sample,
      chrom: requiredText(row, 'chrom', label),
      start,
      end,
      total_cn: totalCn,
      minor_cn: minorCn,
    }
  })
  const assignments = assignmentsValue.map((assignment, index) => {
    if (!isRecord(assignment)) {
      throw new Error(`Invalid simulation data: mutation assignment ${index} is invalid.`)
    }
    const clone = finiteNumber(assignment.clone, `mutation_assignments[${index}].clone`)
    const trueCcf = finiteNumber(assignment.true_ccf, `mutation_assignments[${index}].true_ccf`)
    if (
      assignment.mutation_id !== index ||
      !Number.isSafeInteger(clone) ||
      clone < 0 ||
      clone >= config.clone_count ||
      trueCcf < 0 ||
      trueCcf > 1
    ) {
      throw new Error(`Invalid simulation data: mutation assignment ${index} is out of range.`)
    }
    const mutation = mutations[index]
    if (
      assignment.gene !== mutation.gene ||
      assignment.chrom !== mutation.chrom ||
      assignment.pos !== mutation.pos ||
      Math.abs(trueCcf - cloneFractions[clone]) > 1e-6
    ) {
      throw new Error(`Invalid simulation data: mutation assignment ${index} does not match its row.`)
    }
    return {
      mutation_id: index,
      clone,
      gene: mutation.gene,
      chrom: mutation.chrom,
      pos: mutation.pos,
      true_ccf: trueCcf,
    }
  })

  if (mutations.length !== config.mutation_count || assignments.length !== mutations.length) {
    throw new Error('Invalid simulation data: mutation data does not match the configuration.')
  }

  return {
    config,
    parents: normalizedParents,
    cloneFractions,
    assignments,
    mutations,
    segments,
  }
}

export async function loadDemoSimulation(): Promise<SimulationBundle> {
  const paths = [
    '/sample-bundle/ground_truth.json',
    '/sample-bundle/mutations.tsv',
    '/sample-bundle/copy_number_segments.tsv',
  ]
  const responses = await Promise.all(paths.map((path) => fetch(path)))
  const failedIndex = responses.findIndex((response) => !response.ok)
  if (failedIndex !== -1) {
    throw new Error(`Could not load the bundled simulation file: ${paths[failedIndex]}`)
  }
  const [groundTruth, mutations, segments] = await Promise.all([
    responses[0].json() as Promise<unknown>,
    responses[1].text(),
    responses[2].text(),
  ])
  return parseSimulationBundle(groundTruth, mutations, segments)
}
