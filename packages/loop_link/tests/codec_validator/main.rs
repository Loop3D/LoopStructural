//! Check frames with the real geoviewer codec.
//!
//! Usage: codec-validator FILE...
//! For each file the program reads the bytes, decodes the frame and validates the object.
//! It prints `OK <path>` or `FAIL <path>: <error>`.
//! Exit code: 0 if all files are valid, 1 if one file is not valid, 2 for a usage error.

use geoviewer_protocol::Frame;

fn main() {
    let paths: Vec<String> = std::env::args().skip(1).collect();
    if paths.is_empty() {
        eprintln!("usage: codec-validator FILE...");
        std::process::exit(2);
    }
    let mut failed = false;
    for path in &paths {
        let result = match std::fs::read(path) {
            Err(error) => Err(format!("cannot read file: {error}")),
            Ok(bytes) => Frame::decode(&bytes)
                .map_err(|error| format!("decode: {error}"))
                .and_then(|frame| {
                    frame
                        .validate_object()
                        .map_err(|error| format!("validate_object: {error}"))
                }),
        };
        match result {
            Ok(()) => println!("OK {path}"),
            Err(message) => {
                failed = true;
                println!("FAIL {path}: {message}");
            }
        }
    }
    std::process::exit(i32::from(failed));
}
