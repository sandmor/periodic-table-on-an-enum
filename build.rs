extern crate json;
extern crate memchr;

use std::fs::File;
use std::env;
use std::path::PathBuf;
use std::io::{Read, Write};
use std::u8;

use json::JsonValue;
use memchr::memchr;

#[derive(Debug)]
struct Record {
    atomic_number: u8,
    symbol: String,
    name: String,
    atomic_mass: String,
    cpk: String,
    electron_configuration: String,
    electronegativity: String,
    atomic_radius: String,
    ionization_energy: String,
    electron_affinity: String,
    oxidation_states: Vec<i8>,
    standard_state: String,
    melting_point: String,
    boiling_point: String,
    density: String,
    group_block: String,
    year_discovered: String,
}

fn parse_dhex(i: &[u8]) -> u8 {
    let mut res;
    if i[0] >= b'a' {
        res = i[0] - b'a' + 10;
    }
    else if i[0] >= b'A' {
        res = i[0] - b'A' + 10;
    }
    else {
        res = i[0] - b'0';
    }
    res <<= 4;
    if i[1] >= b'a' {
        res += i[1] - b'a' + 10;
    }
    else if i[1] >= b'A' {
        res += i[1] - b'A' + 10;
    }
    else {
        res += i[1] - b'0';
    }
    res
}

fn main() {
    let mut oxn = 0;
    let mut elements_path = PathBuf::from(env::var("CARGO_MANIFEST_DIR").unwrap());
    let mut out_file = PathBuf::from(env::var("OUT_DIR").unwrap());
    elements_path.push("PubChemElements_all.json");
    out_file.push("data.rs");
    let mut elements_file = File::open(elements_path).unwrap();
    let mut elements_data = String::new();
    elements_file.read_to_string(&mut elements_data).unwrap();
    let elements_data = match json::parse(&elements_data).unwrap() {
        JsonValue::Object(o) => o,
        _ => panic!("Corrupted JSON file")
    };
    let elements_data = match elements_data.get("Table") {
        Some(JsonValue::Object(o)) => o,
        _ => panic!("Corrupted JSON file")
    };
    let columns = match elements_data.get("Columns") {
        Some(JsonValue::Object(o)) => match o.get("Column") {
            Some(JsonValue::Array(arr)) => {
                let mut columns = Vec::with_capacity(arr.len());
                for col in arr {
                    match col {
                        JsonValue::String(col) => {
                            columns.push(col.to_owned());
                        },
                        JsonValue::Short(col) => {
                            columns.push(col.as_str().to_owned());
                        },
                        _ => panic!("Corrupted JSON file")
                    }
                }
                columns
            },
            _ => panic!("Corrupted JSON file")
        },
        _ => panic!("Corrupted JSON file")
    };
    let mut data = Vec::with_capacity(118);
    let rows = match elements_data.get("Row") {
        Some(JsonValue::Array(arr)) => arr,
        _ => panic!("Corrupted JSON file")
    };
    for row in rows {
        let mut record = Record { atomic_number: 0, symbol: String::new(), name: String::new(), atomic_mass: String::new(),
            cpk: String::new(), electron_configuration: String::new(), electronegativity: String::new(), 
            atomic_radius: String::new(), ionization_energy: String::new(), electron_affinity: String::new(), 
            oxidation_states: Vec::new(), standard_state: String::new(), melting_point: String::new(), 
            boiling_point: String::new(), density: String::new(), group_block: String::new(), year_discovered: String::new() };
        let row = match row {
            JsonValue::Object(o) => o,
            _ => panic!("Corrupted JSON file")
        };
        let row = match row.get("Cell") {
            Some(JsonValue::Array(arr)) => {
                let mut rows = Vec::with_capacity(arr.len());
                for row in arr {
                    match row {
                        JsonValue::String(row) => {
                            rows.push(row.to_owned());
                        },
                        JsonValue::Short(row) => {
                            rows.push(row.as_str().to_owned());
                        },
                        _ => panic!("Corrupted JSON file")
                    }
                }
                rows
            },
            _ => panic!("Corrupted JSON file")
        };
        for (i, v) in row.into_iter().enumerate() {
            if v.is_empty() {
                continue;
            }
            match &columns[i][..] {
                "AtomicNumber" => {
                    record.atomic_number = v.parse().unwrap();
                },
                "Symbol" => {
                    record.symbol = v;
                },
                "Name" => {
                    record.name = v;
                },
                "AtomicMass" => {
                    record.atomic_mass = v;
                },
                "CPKHexColor" => {
                    record.cpk = v;
                },
                "ElectronConfiguration" => {
                    record.electron_configuration = v;
                },
                "Electronegativity" => {
                    record.electronegativity = v;
                },
                "AtomicRadius" => {
                    record.atomic_radius = v;
                },
                "IonizationEnergy" => {
                    record.ionization_energy = v;
                },
                "ElectronAffinity" => {
                    record.electron_affinity = v;
                },
                "OxidationStates" => {
                    for state in v.split(',') {
                        let mut state = state.trim().chars();
                        let mut result = 0i8;
                        let mut n = Vec::with_capacity(1);
                        let mut m = 1;
                        match state.next().unwrap() {
                            '-' => {
                                m = -1;
                            },
                            '+' => {},
                            '0' => {},
                            _ => continue
                        }
                        while let Some(c) = state.next() {
                            if !c.is_ascii_digit() {
                                break;
                            }
                            n.push(c as u8 - b'0');
                        }
                        for n in n.into_iter().rev() {
                            result += n as i8 * m;
                            m *= 10;
                        }
                        oxn += 1;
                        record.oxidation_states.push(result);
                    }
                },
                "StandardState" => {
                    record.standard_state = v;
                },
                "MeltingPoint" => {
                    record.melting_point = v;
                },
                "BoilingPoint" => {
                    record.boiling_point = v;
                },
                "Density" => {
                    record.density = v;
                },
                "GroupBlock" => {
                    record.group_block = v;
                },
                "YearDiscovered" => {
                    record.year_discovered = v;
                },
                _ => {}
            }
        }
        data.push(record);
    }
    data.sort_unstable_by_key(|r| r.atomic_number);
    let mut out_file = File::create(out_file).unwrap();
    out_file.write(b"#[derive(Debug, Copy, Clone, Eq, PartialEq, Ord, PartialOrd, Hash)]\npub enum Element {\n").unwrap();
    for record in data.iter() {
        out_file.write(format!("    {},\n", record.name).as_bytes()).unwrap();
    }
    out_file.write(b"}\n").unwrap();
    let mut symbols = Vec::with_capacity(118);
    out_file.write(b"const SYMBOLS: [&str; 118] = [").unwrap();
    let mut first = true;
    for (i, record) in data.iter().enumerate() {
        if first {
            first = false;
        }
        else {
            out_file.write(b", ").unwrap();
        }
        symbols.push((record.symbol.clone(), i));
        out_file.write(format!("\"{}\"", record.symbol).as_bytes()).unwrap();
    }
    out_file.write(b"];\n").unwrap();
    symbols.sort_unstable_by(|a, b| a.0.cmp(&b.0));
    out_file.write(b"const SYMBOLS_SORTED_ALPHABETICALLY: [(&str, u8); 118] = [").unwrap();
    let mut first = true;
    for symbol in symbols {
        if first {
            first = false;
        }
        else {
            out_file.write(b", ").unwrap();
        }
        out_file.write(format!("(\"{}\", {})", symbol.0, symbol.1).as_bytes()).unwrap();
    }
    out_file.write(b"];\n").unwrap();
    let mut names = Vec::with_capacity(118);
    out_file.write(b"const NAMES: [&str; 118] = [").unwrap();
    first = true;
    for (i, record) in data.iter().enumerate() {
        if first {
            first = false;
        }
        else {
            out_file.write(b", ").unwrap();
        }
        names.push((record.name.clone(), i));
        out_file.write(format!("\"{}\"", record.name).as_bytes()).unwrap();
    }
    out_file.write(b"];\n").unwrap();
    names.sort_unstable_by(|a, b| a.0.cmp(&b.0));
    out_file.write(b"const LOWERCASE_NAMES_SORTED_ALPHABETICALLY: [(&str, u8); 118] = [").unwrap();
    let mut first = true;
    for name in names {
        if first {
            first = false;
        }
        else {
            out_file.write(b", ").unwrap();
        }
        out_file.write(format!("(\"{}\", {})", name.0, name.1).as_bytes()).unwrap();
    }
    out_file.write(b"];\n").unwrap();
    out_file.write(b"const ATOMIC_MASSES: [f32; 118] = [").unwrap();
    first = true;
    for record in data.iter() {
        if first {
            first = false;
        }
        else {
            out_file.write(b", ").unwrap();
        }
        out_file.write(record.atomic_mass.as_bytes()).unwrap();
        if let None = memchr(b'.', record.atomic_mass.as_bytes()) {
            out_file.write(b".").unwrap();
        }
    }
    out_file.write(b"];\n").unwrap();
    out_file.write(b"const CPK: [[u8; 3]; 118] = [").unwrap();
    first = true;
    for record in data.iter() {
        if first {
            first = false;
        }
        else {
            out_file.write(b", ").unwrap();
        }
        if record.cpk.is_empty() {
            out_file.write(b"[0, 0, 0]").unwrap();
            continue;
        }
        let hex = record.cpk.as_bytes();
        let r = parse_dhex(&hex[0..]);
        let g = parse_dhex(&hex[2..]);
        let b;
        if hex.len() >= 6 {
            b = parse_dhex(&hex[4..]);
        }
        else {
            b = 0;
        }
        out_file.write(format!("[{}, {}, {}]", r, g, b).as_bytes()).unwrap();
    }
    out_file.write(b"];\n").unwrap();
    out_file.write(b"const ELECTRON_CONFIGURATIONS: [&str; 118] = [").unwrap();
    first = true;
    for record in data.iter() {
        if first {
            first = false;
        }
        else {
            out_file.write(b", ").unwrap();
        }
        out_file.write(format!("\"{}\"", record.electron_configuration).as_bytes()).unwrap();
    }
    out_file.write(b"];\n").unwrap();
    out_file.write(b"const ELECTRONEGATIVITIES: [f32; 118] = [").unwrap();
    first = true;
    for record in data.iter() {
        if first {
            first = false;
        }
        else {
            out_file.write(b", ").unwrap();
        }
        if record.electronegativity.is_empty() {
            out_file.write(b"0.").unwrap();
            continue;
        }
        out_file.write(record.electronegativity.as_bytes()).unwrap();
        if let None = memchr(b'.', record.electronegativity.as_bytes()) {
            out_file.write(b".").unwrap();
        }
    }
    out_file.write(b"];\n").unwrap();
    out_file.write(b"const ATOMIC_RADIUS: [u16; 118] = [").unwrap();
    first = true;
    for record in data.iter() {
        if first {
            first = false;
        }
        else {
            out_file.write(b", ").unwrap();
        }
        if record.atomic_radius.is_empty() {
            out_file.write(b"0").unwrap();
            continue;
        }
        out_file.write(record.atomic_radius.as_bytes()).unwrap();
    }
    out_file.write(b"];\n").unwrap();
    out_file.write(b"const IONIZATION_ENERGIES: [f32; 118] = [").unwrap();
    first = true;
    for record in data.iter() {
        if first {
            first = false;
        }
        else {
            out_file.write(b", ").unwrap();
        }
        if record.ionization_energy.is_empty() {
            out_file.write(b"0.").unwrap();
            continue;
        }
        out_file.write(record.ionization_energy.as_bytes()).unwrap();
        if let None = memchr(b'.', record.ionization_energy.as_bytes()) {
            out_file.write(b".").unwrap();
        }
    }
    out_file.write(b"];\n").unwrap();
    out_file.write(b"const ELECTRON_AFFINITIES: [f32; 118] = [").unwrap();
    first = true;
    for record in data.iter() {
        if first {
            first = false;
        }
        else {
            out_file.write(b", ").unwrap();
        }
        if record.electron_affinity.is_empty() {
            out_file.write(b"0.").unwrap();
            continue;
        }
        out_file.write(record.electron_affinity.as_bytes()).unwrap();
        if let None = memchr(b'.', record.electron_affinity.as_bytes()) {
            out_file.write(b".").unwrap();
        }
    }
    out_file.write(b"];\n").unwrap();
    out_file.write(format!("const OXIDATION_STATES_DATA: [i8; {}] = [", oxn).as_bytes()).unwrap();
    first = true;
    for record in data.iter() {
        for ox in record.oxidation_states.iter() {
            if first {
                first = false;
            }
            else {
                out_file.write(b", ").unwrap();
            }
            out_file.write(format!("{}", ox).as_bytes()).unwrap();
        }
    }
    out_file.write(b"];\n").unwrap();
    out_file.write(b"const OXIDATION_STATES: [(u8, u8); 118] = [").unwrap();
    first = true;
    let mut oi = 0;
    for record in data.iter() {
        if first {
            first = false;
        }
        else {
            out_file.write(b", ").unwrap();
        }
        out_file.write(format!("({}, {})", oi, record.oxidation_states.len()).as_bytes()).unwrap();
        oi += record.oxidation_states.len();
    }
    out_file.write(b"];\n").unwrap();
    out_file.write(b"const STANDARD_STATES: [StateOfMatter; 118] = [").unwrap();
    first = true;
    for record in data.iter() {
        if first {
            first = false;
        }
        else {
            out_file.write(b", ").unwrap();
        }
        out_file.write(match &record.standard_state.to_lowercase()[..] {
            "solid" => b"StateOfMatter::Solid",
            "liquid" => b"StateOfMatter::Liquid",
            "gas" => b"StateOfMatter::Gas",
            "expected to be a solid" => b"StateOfMatter::Solid",
            "expected to be a liquid" => b"StateOfMatter::Liquid",
            "expected to be a gas" => b"StateOfMatter::Gas",
            p@_ => panic!("Unknown state: {}", p)
        }).unwrap();
    }
    out_file.write(b"];\n").unwrap();
    out_file.write(b"const MELTING_POINTS: [f32; 118] = [").unwrap();
    first = true;
    for record in data.iter() {
        if first {
            first = false;
        }
        else {
            out_file.write(b", ").unwrap();
        }
        if record.melting_point.is_empty() {
            out_file.write(b"0.").unwrap();
            continue;
        }
        out_file.write(record.melting_point.as_bytes()).unwrap();
        if let None = memchr(b'.', record.melting_point.as_bytes()) {
            out_file.write(b".").unwrap();
        }
    }
    out_file.write(b"];\n").unwrap();
    out_file.write(b"const BOILING_POINTS: [f32; 118] = [").unwrap();
    first = true;
    for record in data.iter() {
        if first {
            first = false;
        }
        else {
            out_file.write(b", ").unwrap();
        }
        if record.boiling_point.is_empty() {
            out_file.write(b"0.").unwrap();
            continue;
        }
        out_file.write(record.boiling_point.as_bytes()).unwrap();
        if let None = memchr(b'.', record.boiling_point.as_bytes()) {
            out_file.write(b".").unwrap();
        }
    }
    out_file.write(b"];\n").unwrap();
    out_file.write(b"const DENSITIES: [f32; 118] = [").unwrap();
    first = true;
    for record in data.iter() {
        if first {
            first = false;
        }
        else {
            out_file.write(b", ").unwrap();
        }
        if record.density.is_empty() {
            out_file.write(b"0.").unwrap();
            continue;
        }
        out_file.write(record.density.as_bytes()).unwrap();
        if let None = memchr(b'.', record.density.as_bytes()) {
            out_file.write(b".").unwrap();
        }
    }
    out_file.write(b"];\n").unwrap();
    out_file.write(b"const GROUPS: [GroupBlock; 118] = [").unwrap();
    first = true;
    for record in data.iter() {
        if first {
            first = false;
        }
        else {
            out_file.write(b", ").unwrap();
        }
        out_file.write(match &record.group_block.to_lowercase()[..] {
            "halogen" => b"GroupBlock::Halogen",
            "noble gas" => b"GroupBlock::NobleGas",
            "transition metal" => b"GroupBlock::TransitionMetal",
            "post-transition metal" => b"GroupBlock::PostTransitionMetal",
            "alkali metal" => b"GroupBlock::AlkaliMetal",
            "alkaline earth metal" => b"GroupBlock::AlkalineEarthMetal",
            "metalloid" => b"GroupBlock::Metalloid",
            "nonmetal" => b"GroupBlock::NonMetal",
            "actinide" => b"GroupBlock::Actinide",
            "lanthanide" => b"GroupBlock::Lanthanide",
            p@_ => panic!("Unknown group block: {}", p)
        }).unwrap();
    }
    out_file.write(b"];\n").unwrap();
    out_file.write(b"const YEARS_OF_DISCOVERED: [u16; 118] = [").unwrap();
    first = true;
    for record in data.iter() {
        if first {
            first = false;
        }
        else {
            out_file.write(b", ").unwrap();
        }
        if record.year_discovered.is_empty() || record.year_discovered == "Ancient" {
            out_file.write(b"0").unwrap();
            continue;
        }
        out_file.write(record.year_discovered.as_bytes()).unwrap();
    }
    out_file.write(b"];\n").unwrap();
}