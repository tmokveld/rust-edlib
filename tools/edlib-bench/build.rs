use std::{env, path::PathBuf};

fn main() {
    println!("cargo:rerun-if-env-changed=EDLIB_SOURCE");
    let source = PathBuf::from(
        env::var_os("EDLIB_SOURCE").expect("set EDLIB_SOURCE to the upstream Edlib checkout"),
    );
    let edlib = source.join("edlib/src/edlib.cpp");
    let include = source.join("edlib/include");
    println!("cargo:rerun-if-changed={}", edlib.display());
    println!(
        "cargo:rerun-if-changed={}",
        include.join("edlib.h").display()
    );
    println!("cargo:rerun-if-changed=upstream.cpp");
    cc::Build::new()
        .cpp(true)
        .std("c++14")
        .opt_level(3)
        .define("NDEBUG", None)
        .include(include)
        .file(edlib)
        .file("upstream.cpp")
        .compile("edlib_upstream");
}
