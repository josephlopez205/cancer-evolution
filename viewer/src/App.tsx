import { useEffect, useMemo, useState } from 'react'
import './App.css'
import {
  loadDemoSimulation,
  type SimulationBundle,
  type MutationAssignment,
  type MutationObservation,
} from './lib/simulation'

const cloneColors = ['#335f57', '#c8784d', '#6b77a8', '#c4a24e', '#8d6ca8', '#4b8c91']

interface TreeNodePosition {
  x: number
  y: number
}

function getTreePositions(parents: number[]): { positions: TreeNodePosition[]; height: number } {
  const children = parents.map(() => [] as number[])
  parents.forEach((parent, child) => {
    if (child !== 0) children[parent].push(child)
  })

  const positions = parents.map(() => ({ x: 0, y: 0 }))
  let leafIndex = 0
  const place = (node: number, depth: number): number => {
    const nodeChildren = children[node]
    let y: number
    if (nodeChildren.length === 0) {
      y = 52 + leafIndex * 72
      leafIndex += 1
    } else {
      const childYs = nodeChildren.map((child) => place(child, depth + 1))
      y = childYs.reduce((sum, childY) => sum + childY, 0) / childYs.length
    }
    positions[node] = { x: 54 + depth * 138, y }
    return y
  }
  place(0, 0)
  return { positions, height: Math.max(180, 104 + (leafIndex - 1) * 72) }
}

function BrandMark() {
  return (
    <svg aria-hidden="true" className="brand-mark" viewBox="0 0 34 34">
      <path d="M7 17h8m4 0h8M15 9l4 8-4 8" />
      <circle cx="7" cy="17" r="3" />
      <circle cx="27" cy="17" r="3" />
      <circle cx="15" cy="9" r="3" />
      <circle cx="15" cy="25" r="3" />
    </svg>
  )
}

function Icon({ name }: { name: 'sample' | 'cohort' | 'benchmark' }) {
  const paths = {
    sample: <><path d="M5 19.5 12 23l7-3.5V8.5L12 5 5 8.5z" /><path d="m5 8.5 7 3.5 7-3.5M12 12v11" /></>,
    cohort: <><circle cx="9" cy="8" r="3" /><circle cx="17" cy="9" r="2.5" /><path d="M3.5 20c.5-3.2 2.3-5 5.5-5s5 1.8 5.5 5M14 15c3.8-.5 6 1.2 6.5 4" /></>,
    benchmark: <><path d="M5 20V12m7 8V5m7 15v-5" /><path d="M3 21h18" /></>,
  }
  return (
    <svg aria-hidden="true" className="nav-icon" viewBox="0 0 24 24">
      {paths[name]}
    </svg>
  )
}

function Sidebar() {
  return (
    <aside className="sidebar">
      <a className="brand" href="#" aria-label="CloneTrace home">
        <BrandMark />
        <span>CloneTrace</span>
      </a>
      <div className="workspace-label">WORKSPACE</div>
      <nav aria-label="Main navigation">
        <a className="nav-item active" href="#sample" aria-current="page">
          <Icon name="sample" />
          <span>Simulation sample</span>
        </a>
        <span className="nav-item unavailable" aria-disabled="true" title="Available when cohort analysis is implemented">
          <Icon name="cohort" />
          <span>Cohort</span>
          <span className="nav-soon">Later</span>
        </span>
        <span className="nav-item unavailable" aria-disabled="true" title="Available when benchmark runs are implemented">
          <Icon name="benchmark" />
          <span>Benchmarks</span>
          <span className="nav-soon">Later</span>
        </span>
      </nav>
      <div className="sidebar-bottom">
        <div className="sidebar-project">PROJECT</div>
        <div className="project-name"><span className="project-dot" />CloneTrace / BRCA</div>
        <div className="project-note">Research preview · Phase 2</div>
      </div>
    </aside>
  )
}

function StatCard({
  label,
  value,
  detail,
  icon,
  color,
}: {
  label: string
  value: string
  detail: string
  icon: string
  color: string
}) {
  return (
    <article className="stat-card">
      <div className={`stat-icon ${color}`} aria-hidden="true">{icon}</div>
      <div className="stat-copy">
        <span className="stat-label">{label}</span>
        <strong className="stat-value">{value}</strong>
        <span className="stat-detail">{detail}</span>
      </div>
    </article>
  )
}

function CloneTree({ bundle }: { bundle: SimulationBundle }) {
  const { positions, height } = useMemo(() => getTreePositions(bundle.parents), [bundle.parents])
  const counts = bundle.assignments.reduce(
    (result, assignment) => {
      result[assignment.clone] += 1
      return result
    },
    bundle.parents.map(() => 0),
  )

  return (
    <div className="tree-wrap">
      <svg
        className="tree-svg"
        viewBox={`0 0 650 ${height}`}
        role="img"
        aria-label={`Ground-truth clone tree with ${bundle.parents.length} clones`}
      >
        {bundle.parents.map((parent, child) => {
          if (child === 0) return null
          const source = positions[parent]
          const target = positions[child]
          const midpoint = (source.x + target.x) / 2
          return (
            <path
              className="tree-edge"
              key={`edge-${parent}-${child}`}
              d={`M ${source.x + 14} ${source.y} C ${midpoint} ${source.y}, ${midpoint} ${target.y}, ${target.x - 14} ${target.y}`}
            />
          )
        })}
        {positions.map((position, clone) => (
          <g key={`clone-${clone}`}>
            <circle
              className="tree-node"
              cx={position.x}
              cy={position.y}
              r="13"
              fill={cloneColors[clone % cloneColors.length]}
            />
            <text className="tree-node-number" x={position.x} y={position.y + 4}>
              {clone}
            </text>
            <text className="tree-node-label" x={position.x + 23} y={position.y - 5}>
              Clone {clone}
            </text>
            <text className="tree-node-meta" x={position.x + 23} y={position.y + 12}>
              {(bundle.cloneFractions[clone] * 100).toFixed(1)}% · {counts[clone]} variants
            </text>
          </g>
        ))}
      </svg>
    </div>
  )
}

function CopyNumberPanel({ bundle }: { bundle: SimulationBundle }) {
  if (bundle.segments.length === 0) {
    return <p className="empty-note">No copy-number events were generated for this sample.</p>
  }
  return (
    <div className="segment-list">
      {bundle.segments.map((segment, index) => (
        <div className="segment-row" key={`${segment.chrom}-${segment.start}-${index}`}>
          <span className="chromosome-tag">chr{segment.chrom}</span>
          <span className="segment-range">
            {segment.start.toLocaleString()}–{segment.end.toLocaleString()}
          </span>
          <span className="cn-value">
            <b>{segment.total_cn}</b>
            <small>total CN</small>
          </span>
          <span className="cn-value minor">
            <b>{segment.minor_cn}</b>
            <small>minor</small>
          </span>
        </div>
      ))}
    </div>
  )
}

function MutationTable({
  mutations,
  assignments,
  selectedIndex,
  onSelect,
}: {
  mutations: MutationObservation[]
  assignments: MutationAssignment[]
  selectedIndex: number
  onSelect: (index: number) => void
}) {
  return (
    <div className="table-scroll">
      <table className="mutation-table">
        <thead>
          <tr>
            <th>Mutation</th>
            <th>Clone</th>
            <th>Read counts</th>
            <th>VAF</th>
            <th>True CCF</th>
          </tr>
        </thead>
        <tbody>
          {mutations.map((mutation, index) => {
            const assignment = assignments[index]
            const totalReads = mutation.ref_count + mutation.alt_count
            const vaf = totalReads ? mutation.alt_count / totalReads : 0
            return (
              <tr
                aria-selected={selectedIndex === index}
                className={selectedIndex === index ? 'selected-row' : ''}
                key={`${mutation.chrom}-${mutation.pos}-${index}`}
                onClick={() => onSelect(index)}
                onKeyDown={(event) => {
                  if (event.key === 'Enter' || event.key === ' ') onSelect(index)
                }}
                tabIndex={0}
              >
                <td>
                  <span className="gene-name">{mutation.gene}</span>
                  <span className="variant-name">
                    chr{mutation.chrom}:{mutation.pos} {mutation.ref}→{mutation.alt}
                  </span>
                </td>
                <td>
                  <span className="clone-pill">
                    <i style={{ backgroundColor: cloneColors[assignment.clone % cloneColors.length] }} />
                    C{assignment.clone}
                  </span>
                </td>
                <td className="read-count">{mutation.alt_count}<span> alt</span> / {mutation.ref_count}<span> ref</span></td>
                <td>{(vaf * 100).toFixed(1)}%</td>
                <td>{(assignment.true_ccf * 100).toFixed(1)}%</td>
              </tr>
            )
          })}
        </tbody>
      </table>
    </div>
  )
}

function SelectedMutation({
  mutation,
  assignment,
}: {
  mutation: MutationObservation
  assignment: MutationAssignment
}) {
  const total = mutation.ref_count + mutation.alt_count
  const alternatePercent = total === 0 ? 0 : (mutation.alt_count / total) * 100
  return (
    <div className="selected-detail">
      <div className="selected-gene">{mutation.gene}</div>
      <div className="selected-position">
        chr{mutation.chrom}:{mutation.pos} <span>{mutation.ref}→{mutation.alt}</span>
      </div>
      <div className="read-legend">
        <span><i className="legend-ref" /> Reference reads <b>{mutation.ref_count}</b></span>
        <span><i className="legend-alt" /> Alternate reads <b>{mutation.alt_count}</b></span>
      </div>
      <div
        className="read-bar"
        role="img"
        aria-label={`${mutation.alt_count} alternate reads and ${mutation.ref_count} reference reads`}
      >
        <span className="read-bar-alt" style={{ width: `${alternatePercent}%` }} />
        <span className="read-bar-ref" style={{ width: `${100 - alternatePercent}%` }} />
      </div>
      <div className="selected-stats">
        <div><span>Observed VAF</span><b>{(alternatePercent).toFixed(1)}%</b></div>
        <div><span>True cancer cell fraction</span><b>{(assignment.true_ccf * 100).toFixed(1)}%</b></div>
        <div><span>Assigned clone</span><b>C{assignment.clone}</b></div>
      </div>
      <p className="detail-disclaimer">True CCF and clone assignment are simulator ground truth, not an inference.</p>
    </div>
  )
}

function App() {
  const [bundle, setBundle] = useState<SimulationBundle | null>(null)
  const [error, setError] = useState<string | null>(null)
  const [selectedIndex, setSelectedIndex] = useState(0)

  useEffect(() => {
    let mounted = true
    loadDemoSimulation()
      .then((result) => {
        if (mounted) setBundle(result)
      })
      .catch((reason: unknown) => {
        if (mounted) {
          setError(reason instanceof Error ? reason.message : 'The simulation sample could not be loaded.')
        }
      })
    return () => {
      mounted = false
    }
  }, [])

  if (error) {
    return (
      <main className="app-shell">
        <Sidebar />
        <section className="page-content state-page">
          <div className="error-card" role="alert">
            <span className="error-symbol">!</span>
            <h1>Could not load sample</h1>
            <p>{error}</p>
            <button type="button" className="secondary-button" onClick={() => window.location.reload()}>
              Try again
            </button>
          </div>
        </section>
      </main>
    )
  }
  if (!bundle) {
    return (
      <main className="app-shell">
        <Sidebar />
        <section className="page-content state-page" aria-live="polite">
          <div className="loading-card"><span className="spinner" />Loading simulated sample…</div>
        </section>
      </main>
    )
  }

  const selectedMutation = bundle.mutations[selectedIndex]
  const selectedAssignment = bundle.assignments[selectedIndex]

  return (
    <main className="app-shell">
      <Sidebar />
      <section className="page-content" id="sample">
        <header className="topbar">
          <div className="breadcrumb"><span>DATA EXPLORER</span><b>/</b> SIMULATION</div>
          <div className="topbar-right">
            <span className="data-state"><i /> Local demo data</span>
            <span className="avatar" aria-label="CloneTrace">CT</span>
          </div>
        </header>

        <div className="page-main">
          <div className="page-heading">
            <div>
              <div className="eyebrow">PHASE 02 <span>·</span> SIMULATION BENCHMARK</div>
              <h1>Simulation sample</h1>
              <p>Explore a generated clonal architecture and its observed sequencing reads.</p>
            </div>
            <div className="sample-actions">
              <span className="truth-badge"><i /> Simulated ground truth</span>
              <button type="button" className="secondary-button" disabled title="Export will be available with the bundle contract">
                Export sample
              </button>
            </div>
          </div>

          <div className="sample-meta">
            <span><b>Sample</b> {bundle.config.sample}</span>
            <span><b>Seed</b> {bundle.config.seed}</span>
            <span><b>Topology</b> {bundle.config.topology}</span>
            <span><b>Copy-number events</b> {bundle.segments.length}</span>
          </div>

          <section className="stat-grid" aria-label="Simulation summary">
            <StatCard label="Clones" value={String(bundle.config.clone_count)} detail="in ground-truth tree" icon="◉" color="green" />
            <StatCard label="Mutations" value={String(bundle.mutations.length)} detail="simulated variants" icon="⌁" color="orange" />
            <StatCard label="Read depth" value={`${bundle.config.depth}×`} detail="target coverage" icon="▥" color="blue" />
            <StatCard label="Tumor purity" value={`${Math.round(bundle.config.purity * 100)}%`} detail="configured purity" icon="◒" color="gold" />
          </section>

          <section className="content-grid">
            <article className="panel tree-panel">
              <div className="panel-heading">
                <div>
                  <span className="panel-kicker">CLONAL ARCHITECTURE</span>
                  <h2>Ground-truth clone tree</h2>
                </div>
                <span className="panel-count">{bundle.config.clone_count} clones</span>
              </div>
              <CloneTree bundle={bundle} />
              <div className="tree-footnote">
                <span><i className="tree-legend-dot" /> Clone prevalence shown as fraction of tumor cells</span>
                <span>Root → descendants</span>
              </div>
            </article>

            <article className="panel selected-panel">
              <div className="panel-heading">
                <div>
                  <span className="panel-kicker">SELECTED VARIANT</span>
                  <h2>Observed reads</h2>
                </div>
                <span className="read-depth-tag">{selectedMutation.ref_count + selectedMutation.alt_count} reads</span>
              </div>
              <SelectedMutation mutation={selectedMutation} assignment={selectedAssignment} />
            </article>

            <article className="panel mutation-panel">
              <div className="panel-heading">
                <div>
                  <span className="panel-kicker">VARIANT OBSERVATIONS</span>
                  <h2>Mutation read counts</h2>
                  <p className="panel-description">Select a row to inspect alternate and reference reads.</p>
                </div>
                <span className="panel-count">{bundle.mutations.length} variants</span>
              </div>
              <MutationTable
                mutations={bundle.mutations}
                assignments={bundle.assignments}
                selectedIndex={selectedIndex}
                onSelect={setSelectedIndex}
              />
            </article>

            <article className="panel cn-panel">
              <div className="panel-heading">
                <div>
                  <span className="panel-kicker">COPY NUMBER</span>
                  <h2>Simulated events</h2>
                </div>
                <span className="panel-count">{bundle.segments.length} events</span>
              </div>
              <CopyNumberPanel bundle={bundle} />
            </article>
          </section>

          <footer className="page-footer">
            <span><b>CloneTrace preview</b> · This sample was generated by the simulator.</span>
            <span>Not a clinical or diagnostic result.</span>
          </footer>
        </div>
      </section>
    </main>
  )
}

export default App
