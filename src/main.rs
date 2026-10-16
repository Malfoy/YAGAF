use clap::Parser;
fn main() -> anyhow::Result<()> {
    yagaf::pipeline::run(&yagaf::pipeline::Args::parse())
}
