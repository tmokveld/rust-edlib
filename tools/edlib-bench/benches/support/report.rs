use super::{Case, Outcome};
use serde::Deserialize;
use std::{
    env, fmt,
    fs::{self, File},
    io::{self, BufWriter, Write},
    path::{Path, PathBuf},
};

pub(super) struct Table {
    home: PathBuf,
    baseline: Option<Baseline>,
    out: BufWriter<File>,
    started: bool,
}

struct Baseline {
    home: PathBuf,
    revision: String,
}

impl Table {
    pub(super) fn from_env() -> io::Result<Option<Self>> {
        let Some(summary) = env::var_os("EDLIB_BENCH_SUMMARY") else {
            return Ok(None);
        };
        let home = env::var_os("CRITERION_HOME").expect("run.sh sets CRITERION_HOME");
        let baseline = env::var_os("EDLIB_BENCH_BASELINE_HOME").map(|home| Baseline {
            home: home.into(),
            revision: env::var("EDLIB_BENCH_BASELINE_REV").expect("run.sh sets the baseline revision"),
        });
        Ok(Some(Self {
            home: home.into(),
            baseline,
            out: BufWriter::new(File::create(summary)?),
            started: false,
        }))
    }

    pub(super) fn case_dir(&self, group: &str, case: &str) -> PathBuf {
        case_dir(&self.home, group, case)
    }

    pub(super) fn write_row(
        &mut self,
        group: &str,
        case: &Case<'_>,
        outcome: &Outcome,
        rust: Option<f64>,
        edlib: Option<f64>,
    ) -> io::Result<()> {
        let out = &mut self.out;
        if !self.started {
            writeln!(out, "\nCriterion median ns/call. Ratios above 1× favor current Rust.")?;
            if let Some(base) = &self.baseline {
                writeln!(out, "\nRust baseline: `{}`.", base.revision)?;
            }
            write!(out, "\n| Group | Input | Search | Task | Query bytes | Target bytes | Bound | Score | Rust ns/call | Edlib ns/call | Edlib / Rust |")?;
            if self.baseline.is_some() {
                write!(out, " Baseline ns/call | Baseline / Rust |")?;
            }
            write!(out, "\n|---|---|---|---|---:|---:|---:|---:|---:|---:|---:|")?;
            if self.baseline.is_some() {
                write!(out, "---:|---:|")?;
            }
            writeln!(out)?;
            self.started = true;
        }
        write!(
            out,
            "| {group} | {} | {} | {} | {} | {} | ",
            case.input, case.mode, case.task, case.query_len, case.target_len,
        )?;
        match outcome.cutoff {
            Some(k) => write!(out, "{k} | ")?,
            None => write!(out, "unbounded | ")?,
        }
        match outcome.score {
            Some(score) => write!(out, "{score} | ")?,
            None => write!(out, "no match | ")?,
        }
        write!(out, "{} | {} | {} |", Timing(rust), Timing(edlib), Ratio(edlib, rust))?;
        if let Some(base) = &self.baseline {
            // Both runs use identical Rust IDs and filters. No Rust timing
            // here means the case was filtered out of the baseline as well.
            let baseline = rust.map(|_| {
                let dir = case_dir(&base.home, group, &case.to_string());
                median(&estimate_path(&dir, "rust"))
            });
            write!(out, " {} | {} |", Timing(baseline), Ratio(baseline, rust))?;
        }
        writeln!(out)
    }

    pub(super) fn finish(mut self) -> io::Result<()> {
        self.out.flush()
    }
}

fn case_dir(home: &Path, group: &str, case: &str) -> PathBuf {
    // Criterion 0.8.2 replaces slashes in each ID component with underscores.
    // Our components are unique, short, and otherwise filename-safe.
    home.join(group).join(case.replace('/', "_"))
}

pub(super) fn estimate_path(dir: &Path, engine: &str) -> PathBuf {
    dir.join(engine).join("new/estimates.json")
}

struct Timing(Option<f64>);

impl fmt::Display for Timing {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self.0 {
            Some(ns) => write!(f, "{ns:.0}"),
            None => f.write_str("—"),
        }
    }
}

struct Ratio(Option<f64>, Option<f64>);

impl fmt::Display for Ratio {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self.0.zip(self.1) {
            Some((other, rust)) => write!(f, "{:.2}×", other / rust),
            None => f.write_str("—"),
        }
    }
}

pub(super) fn clear_estimate(path: &Path) {
    if let Err(err) = fs::remove_file(path) {
        if err.kind() != io::ErrorKind::NotFound {
            panic!("removing old estimate {}: {err}", path.display());
        }
    }
}

pub(super) fn median(path: &Path) -> f64 {
    #[derive(Deserialize)]
    struct Estimates {
        median: Estimate,
    }
    #[derive(Deserialize)]
    struct Estimate {
        point_estimate: f64,
    }

    let estimates: Estimates = fs::read(path)
        .and_then(|bytes| Ok(serde_json::from_slice(&bytes)?))
        .unwrap_or_else(|err| panic!("reading new estimate {}: {err}", path.display()));
    let ns = estimates.median.point_estimate;
    assert!(ns.is_finite() && ns > 0.0, "invalid median in {}", path.display());
    ns
}
