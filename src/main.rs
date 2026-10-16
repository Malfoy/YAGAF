use clap::Parser;
fn main() -> anyhow::Result<()> {
    gapfil::pipeline::run(&gapfil::pipeline::Args::parse())
}
