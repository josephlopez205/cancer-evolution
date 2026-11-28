import { describe, expect, it } from 'vitest'
import { parseSimulationBundle } from './simulation'

const groundTruth = {
  format_version: '1',
  config: {
    seed: 1,
    depth: 30,
    purity: 0.5,
    clone_count: 2,
    topology: 'linear',
    mutation_count: 1,
    cn_event_probability: 0.5,
    sample: 'SIM-1',
  },
  clone_tree: {
    parents: [0, 0],
    clone_fractions: [0.6, 0.4],
  },
  mutation_assignments: [
    {
      mutation_id: 0,
      clone: 1,
      gene: 'TP53',
      chrom: '17',
      pos: 100,
      true_ccf: 0.4,
    },
  ],
}
const mutations =
  'sample\tgene\tchrom\tpos\tref\talt\tref_count\talt_count\n' +
  'SIM-1\tTP53\t17\t100\tA\tC\t25\t5\n'
const segments = 'sample\tchrom\tstart\tend\ttotal_cn\tminor_cn\n'

describe('parseSimulationBundle', () => {
  it('parses valid simulator output', () => {
    const bundle = parseSimulationBundle(groundTruth, mutations, segments)

    expect(bundle.config.sample).toBe('SIM-1')
    expect(bundle.parents).toEqual([0, 0])
    expect(bundle.mutations[0].alt_count).toBe(5)
    expect(bundle.assignments[0].clone).toBe(1)
    expect(bundle.segments).toEqual([])
  })

  it('rejects unsupported bundle versions', () => {
    expect(() =>
      parseSimulationBundle({ ...groundTruth, format_version: '2' }, mutations, segments),
    ).toThrow('Unsupported or invalid simulation bundle version.')
  })

  it('rejects mismatched mutation counts', () => {
    expect(() =>
      parseSimulationBundle(
        {
          ...groundTruth,
          config: { ...groundTruth.config, mutation_count: 2 },
        },
        mutations,
        segments,
      ),
    ).toThrow('mutation data does not match the configuration')
  })

  it('rejects assignments that disagree with observed mutation rows', () => {
    const invalidGroundTruth = {
      ...groundTruth,
      mutation_assignments: [
        {
          ...groundTruth.mutation_assignments[0],
          gene: 'PIK3CA',
        },
      ],
    }

    expect(() => parseSimulationBundle(invalidGroundTruth, mutations, segments)).toThrow(
      'does not match its row',
    )
  })
})
