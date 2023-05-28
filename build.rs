extern crate json;

use std::collections::HashMap;
use std::env;
use std::fs::File;
use std::io::{Read, Write};
use std::path::PathBuf;
use std::u8;

use json::JsonValue;

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
    } else if i[0] >= b'A' {
        res = i[0] - b'A' + 10;
    } else {
        res = i[0] - b'0';
    }
    res <<= 4;
    if i[1] >= b'a' {
        res += i[1] - b'a' + 10;
    } else if i[1] >= b'A' {
        res += i[1] - b'A' + 10;
    } else {
        res += i[1] - b'0';
    }
    res
}

#[derive(Copy, Clone)]
pub struct ElectronicConfiguration {
    s: [u8; 7],
    p: [u8; 6],
    d: [u8; 4],
    f: [u8; 2],
}

#[derive(Clone)]
enum EC {
    Unparsed(String),
    Parsed(ElectronicConfiguration),
}

fn parse_number(string: &[u8]) -> (usize, &[u8]) {
    let mut digits = Vec::new();
    let mut end = 0;
    while string.len() > end && (string[end] >= b'0' && string[end] <= b'9') {
        digits.push(string[end] - b'0');
        end += 1;
    }
    let mut result = 0;
    let mut m = 1;
    for d in digits.into_iter().rev() {
        result += d as usize * m;
        m *= 10;
    }
    (result, &string[end..])
}

fn get_ec(symbol: &str, ec: &mut HashMap<String, EC>) -> ElectronicConfiguration {
    match ec.get(symbol).unwrap().clone() {
        EC::Unparsed(s) => {
            eprintln!("{}", s);
            let mut string = s.as_bytes();
            let mut result;
            if string[0] == b'[' {
                string = &string[1..];
                if string[1] == b']' {
                    result = get_ec(std::str::from_utf8(&string[..1]).unwrap(), ec);
                    string = &string[1..];
                } else {
                    result = get_ec(std::str::from_utf8(&string[..2]).unwrap(), ec);
                    string = &string[2..];
                }
                string = &string[1..];
            } else {
                result = ElectronicConfiguration {
                    s: [0; 7],
                    p: [0; 6],
                    d: [0; 4],
                    f: [0; 2],
                };
            }
            while !string.is_empty() {
                if string[0] == b' ' {
                    string = &string[1..];
                }
                if string[0] == b'(' {
                    // Message
                    break;
                }
                let (index, s) = parse_number(string);
                string = s;
                let t = string[0];
                string = &string[1..];
                let (value, s) = parse_number(string);
                string = s;
                match t {
                    b's' => {
                        result.s[index - 1] = value as u8;
                    }
                    b'p' => {
                        result.p[index - 2] = value as u8;
                    }
                    b'd' => {
                        result.d[index - 3] = value as u8;
                    }
                    b'f' => {
                        result.f[index - 4] = value as u8;
                    }
                    _ => unreachable!(),
                }
            }
            ec.insert(symbol.to_owned(), EC::Parsed(result));
            result
        }
        EC::Parsed(c) => c,
    }
}

fn main() {
    let mut ec = HashMap::new();
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
        _ => panic!("Corrupted JSON file"),
    };
    let elements_data = match elements_data.get("Table") {
        Some(JsonValue::Object(o)) => o,
        _ => panic!("Corrupted JSON file"),
    };
    let columns = match elements_data.get("Columns") {
        Some(JsonValue::Object(o)) => match o.get("Column") {
            Some(JsonValue::Array(arr)) => {
                let mut columns = Vec::with_capacity(arr.len());
                for col in arr {
                    match col {
                        JsonValue::String(col) => {
                            columns.push(col.to_owned());
                        }
                        JsonValue::Short(col) => {
                            columns.push(col.as_str().to_owned());
                        }
                        _ => panic!("Corrupted JSON file"),
                    }
                }
                columns
            }
            _ => panic!("Corrupted JSON file"),
        },
        _ => panic!("Corrupted JSON file"),
    };
    let mut data = Vec::with_capacity(118);
    let rows = match elements_data.get("Row") {
        Some(JsonValue::Array(arr)) => arr,
        _ => panic!("Corrupted JSON file"),
    };
    for row in rows {
        let mut record = Record {
            atomic_number: 0,
            symbol: String::new(),
            name: String::new(),
            atomic_mass: String::new(),
            cpk: String::new(),
            electron_configuration: String::new(),
            electronegativity: String::new(),
            atomic_radius: String::new(),
            ionization_energy: String::new(),
            electron_affinity: String::new(),
            oxidation_states: Vec::new(),
            standard_state: String::new(),
            melting_point: String::new(),
            boiling_point: String::new(),
            density: String::new(),
            group_block: String::new(),
            year_discovered: String::new(),
        };
        let row = match row {
            JsonValue::Object(o) => o,
            _ => panic!("Corrupted JSON file"),
        };
        let row = match row.get("Cell") {
            Some(JsonValue::Array(arr)) => {
                let mut rows = Vec::with_capacity(arr.len());
                for row in arr {
                    match row {
                        JsonValue::String(row) => {
                            rows.push(row.to_owned());
                        }
                        JsonValue::Short(row) => {
                            rows.push(row.as_str().to_owned());
                        }
                        _ => panic!("Corrupted JSON file"),
                    }
                }
                rows
            }
            _ => panic!("Corrupted JSON file"),
        };
        for (i, v) in row.into_iter().enumerate() {
            if v.is_empty() {
                continue;
            }
            match &columns[i][..] {
                "AtomicNumber" => {
                    record.atomic_number = v.parse().unwrap();
                }
                "Symbol" => {
                    record.symbol = v;
                }
                "Name" => {
                    record.name = v;
                }
                "AtomicMass" => {
                    record.atomic_mass = v;
                }
                "CPKHexColor" => {
                    record.cpk = v;
                }
                "ElectronConfiguration" => {
                    record.electron_configuration = v;
                }
                "Electronegativity" => {
                    record.electronegativity = v;
                }
                "AtomicRadius" => {
                    record.atomic_radius = v;
                }
                "IonizationEnergy" => {
                    record.ionization_energy = v;
                }
                "ElectronAffinity" => {
                    record.electron_affinity = v;
                }
                "OxidationStates" => {
                    for state in v.split(',') {
                        let mut state = state.trim().chars();
                        let mut result = 0i8;
                        let mut n = Vec::with_capacity(1);
                        let mut m = 1;
                        match state.next().unwrap() {
                            '-' => {
                                m = -1;
                            }
                            '+' => {}
                            '0' => {}
                            _ => continue,
                        }
                        for c in state {
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
                }
                "StandardState" => {
                    record.standard_state = v;
                }
                "MeltingPoint" => {
                    record.melting_point = v;
                }
                "BoilingPoint" => {
                    record.boiling_point = v;
                }
                "Density" => {
                    record.density = v;
                }
                "GroupBlock" => {
                    record.group_block = v;
                }
                "YearDiscovered" => {
                    record.year_discovered = v;
                }
                _ => {}
            }
        }
        ec.insert(
            record.symbol.clone(),
            EC::Unparsed(record.electron_configuration.to_owned()),
        );
        data.push(record);
    }
    data.sort_unstable_by_key(|r| r.atomic_number);
    let mut out_file = File::create(out_file).unwrap();
    out_file.write_all(b"#[derive(Debug, Copy, Clone, Eq, PartialEq, Ord, PartialOrd, Hash)]\n#[non_exhaustive]\npub enum Element {\n").unwrap();
    for record in data.iter() {
        out_file
            .write_all(format!("    {},\n", record.name).as_bytes())
            .unwrap();
    }
    out_file.write_all(b"}\n").unwrap();
    let mut symbols = Vec::with_capacity(118);
    out_file
        .write_all(b"const SYMBOLS: [&str; 118] = [")
        .unwrap();
    let mut first = true;
    for (i, record) in data.iter().enumerate() {
        if first {
            first = false;
        } else {
            out_file.write_all(b", ").unwrap();
        }
        symbols.push((record.symbol.clone(), i));
        out_file
            .write_all(format!("\"{}\"", record.symbol).as_bytes())
            .unwrap();
    }
    out_file.write_all(b"];\n").unwrap();
    symbols.sort_unstable_by(|a, b| a.0.cmp(&b.0));
    out_file
        .write_all(b"const SYMBOLS_SORTED_ALPHABETICALLY: [(&str, u8); 118] = [")
        .unwrap();
    let mut first = true;
    for symbol in symbols {
        if first {
            first = false;
        } else {
            out_file.write_all(b", ").unwrap();
        }
        out_file
            .write_all(format!("(\"{}\", {})", symbol.0, symbol.1).as_bytes())
            .unwrap();
    }
    out_file.write_all(b"];\n").unwrap();
    let mut names = Vec::with_capacity(118);
    out_file.write_all(b"const NAMES: [&str; 118] = [").unwrap();
    first = true;
    for (i, record) in data.iter().enumerate() {
        if first {
            first = false;
        } else {
            out_file.write_all(b", ").unwrap();
        }
        names.push((record.name.to_lowercase(), i));
        out_file
            .write_all(format!("\"{}\"", record.name).as_bytes())
            .unwrap();
    }
    out_file.write_all(b"];\n").unwrap();
    names.sort_unstable_by(|a, b| a.0.cmp(&b.0));
    out_file
        .write_all(b"const LOWERCASE_NAMES_SORTED_ALPHABETICALLY: [(&str, u8); 118] = [")
        .unwrap();
    let mut first = true;
    for name in names {
        if first {
            first = false;
        } else {
            out_file.write_all(b", ").unwrap();
        }
        out_file
            .write_all(format!("(\"{}\", {})", name.0, name.1).as_bytes())
            .unwrap();
    }
    out_file.write_all(b"];\n").unwrap();
    out_file
        .write_all(b"const ATOMIC_MASSES: [f64; 118] = [")
        .unwrap();
    first = true;
    for record in data.iter() {
        if first {
            first = false;
        } else {
            out_file.write_all(b", ").unwrap();
        }
        out_file.write_all(record.atomic_mass.as_bytes()).unwrap();
        if record.atomic_mass.find('.').is_none() {
            out_file.write_all(b".").unwrap();
        }
    }
    out_file.write_all(b"];\n").unwrap();
    out_file
        .write_all(b"const CPK: [[u8; 3]; 118] = [")
        .unwrap();
    first = true;
    for record in data.iter() {
        if first {
            first = false;
        } else {
            out_file.write_all(b", ").unwrap();
        }
        if record.cpk.is_empty() {
            out_file.write_all(b"[0, 0, 0]").unwrap();
            continue;
        }
        let hex = record.cpk.as_bytes();
        let r = parse_dhex(&hex[0..]);
        let g = parse_dhex(&hex[2..]);
        let b = if hex.len() >= 6 {
            parse_dhex(&hex[4..])
        } else {
            0
        };
        out_file
            .write_all(format!("[{}, {}, {}]", r, g, b).as_bytes())
            .unwrap();
    }
    out_file.write_all(b"];\n").unwrap();
    out_file
        .write_all(b"const ELECTRON_CONFIGURATIONS: [&str; 118] = [")
        .unwrap();
    first = true;
    for record in data.iter() {
        if first {
            first = false;
        } else {
            out_file.write_all(b", ").unwrap();
        }
        out_file
            .write_all(format!("\"{}\"", record.electron_configuration).as_bytes())
            .unwrap();
    }
    out_file.write_all(b"];\n").unwrap();
    out_file
        .write_all(b"const ELECTRONIC_CONFIGURATION_PARSED: [ElectronicConfiguration; 118] = [")
        .unwrap();
    first = true;
    for record in data.iter() {
        if first {
            first = false;
        } else {
            out_file.write_all(b", ").unwrap();
        }
        let e = get_ec(&record.symbol, &mut ec);
        out_file
            .write_all(format!("ElectronicConfiguration {} s: [", "{").as_bytes())
            .unwrap();
        let mut inner_first = true;
        for i in e.s.iter() {
            if inner_first {
                inner_first = false;
            } else {
                out_file.write_all(b", ").unwrap();
            }
            out_file.write_all(format!("{}", i).as_bytes()).unwrap();
        }
        out_file.write_all(b"], p: [").unwrap();
        inner_first = true;
        for i in e.p.iter() {
            if inner_first {
                inner_first = false;
            } else {
                out_file.write_all(b", ").unwrap();
            }
            out_file.write_all(format!("{}", i).as_bytes()).unwrap();
        }
        out_file.write_all(b"], d: [").unwrap();
        inner_first = true;
        for i in e.d.iter() {
            if inner_first {
                inner_first = false;
            } else {
                out_file.write_all(b", ").unwrap();
            }
            out_file.write_all(format!("{}", i).as_bytes()).unwrap();
        }
        out_file.write_all(b"], f: [").unwrap();
        inner_first = true;
        for i in e.f.iter() {
            if inner_first {
                inner_first = false;
            } else {
                out_file.write_all(b", ").unwrap();
            }
            out_file.write_all(format!("{}", i).as_bytes()).unwrap();
        }
        out_file.write_all(b"] }").unwrap();
    }
    out_file.write_all(b"];\n").unwrap();
    out_file
        .write_all(b"const ELECTRONEGATIVITIES: [f32; 118] = [")
        .unwrap();
    first = true;
    for record in data.iter() {
        if first {
            first = false;
        } else {
            out_file.write_all(b", ").unwrap();
        }
        if record.electronegativity.is_empty() {
            out_file.write_all(b"0.").unwrap();
            continue;
        }
        out_file
            .write_all(record.electronegativity.as_bytes())
            .unwrap();
        if record.electronegativity.find('.').is_none() {
            out_file.write_all(b".").unwrap();
        }
    }
    out_file.write_all(b"];\n").unwrap();
    out_file
        .write_all(b"const ATOMIC_RADIUS: [u16; 118] = [")
        .unwrap();
    first = true;
    for record in data.iter() {
        if first {
            first = false;
        } else {
            out_file.write_all(b", ").unwrap();
        }
        if record.atomic_radius.is_empty() {
            out_file.write_all(b"0").unwrap();
            continue;
        }
        out_file.write_all(record.atomic_radius.as_bytes()).unwrap();
    }
    out_file.write_all(b"];\n").unwrap();
    out_file
        .write_all(b"const IONIZATION_ENERGIES: [f32; 118] = [")
        .unwrap();
    first = true;
    for record in data.iter() {
        if first {
            first = false;
        } else {
            out_file.write_all(b", ").unwrap();
        }
        if record.ionization_energy.is_empty() {
            out_file.write_all(b"0.").unwrap();
            continue;
        }
        out_file
            .write_all(record.ionization_energy.as_bytes())
            .unwrap();
        if record.ionization_energy.find('.').is_none() {
            out_file.write_all(b".").unwrap();
        }
    }
    out_file.write_all(b"];\n").unwrap();
    out_file
        .write_all(b"const ELECTRON_AFFINITIES: [f32; 118] = [")
        .unwrap();
    first = true;
    for record in data.iter() {
        if first {
            first = false;
        } else {
            out_file.write_all(b", ").unwrap();
        }
        if record.electron_affinity.is_empty() {
            out_file.write_all(b"0.").unwrap();
            continue;
        }
        out_file
            .write_all(record.electron_affinity.as_bytes())
            .unwrap();
        if record.electron_affinity.find('.').is_none() {
            out_file.write_all(b".").unwrap();
        }
    }
    out_file.write_all(b"];\n").unwrap();
    out_file
        .write_all(format!("const OXIDATION_STATES_DATA: [i8; {}] = [", oxn).as_bytes())
        .unwrap();
    first = true;
    for record in data.iter() {
        for ox in record.oxidation_states.iter() {
            if first {
                first = false;
            } else {
                out_file.write_all(b", ").unwrap();
            }
            out_file.write_all(format!("{}", ox).as_bytes()).unwrap();
        }
    }
    out_file.write_all(b"];\n").unwrap();
    out_file
        .write_all(b"const OXIDATION_STATES: [(u8, u8); 118] = [")
        .unwrap();
    first = true;
    let mut oi = 0;
    for record in data.iter() {
        if first {
            first = false;
        } else {
            out_file.write_all(b", ").unwrap();
        }
        out_file
            .write_all(format!("({}, {})", oi, record.oxidation_states.len()).as_bytes())
            .unwrap();
        oi += record.oxidation_states.len();
    }
    out_file.write_all(b"];\n").unwrap();
    out_file
        .write_all(b"const STANDARD_STATES: [StateOfMatter; 118] = [")
        .unwrap();
    first = true;
    for record in data.iter() {
        if first {
            first = false;
        } else {
            out_file.write_all(b", ").unwrap();
        }
        out_file
            .write_all(match &record.standard_state.to_lowercase()[..] {
                "solid" => b"StateOfMatter::Solid",
                "liquid" => b"StateOfMatter::Liquid",
                "gas" => b"StateOfMatter::Gas",
                "expected to be a solid" => b"StateOfMatter::Solid",
                "expected to be a liquid" => b"StateOfMatter::Liquid",
                "expected to be a gas" => b"StateOfMatter::Gas",
                p => panic!("Unknown state: {}", p),
            })
            .unwrap();
    }
    out_file.write_all(b"];\n").unwrap();
    out_file
        .write_all(b"const MELTING_POINTS: [f32; 118] = [")
        .unwrap();
    first = true;
    for record in data.iter() {
        if first {
            first = false;
        } else {
            out_file.write_all(b", ").unwrap();
        }
        if record.melting_point.is_empty() {
            out_file.write_all(b"0.").unwrap();
            continue;
        }
        out_file.write_all(record.melting_point.as_bytes()).unwrap();
        if record.melting_point.find('.').is_none() {
            out_file.write_all(b".").unwrap();
        }
    }
    out_file.write_all(b"];\n").unwrap();
    out_file
        .write_all(b"const BOILING_POINTS: [f32; 118] = [")
        .unwrap();
    first = true;
    for record in data.iter() {
        if first {
            first = false;
        } else {
            out_file.write_all(b", ").unwrap();
        }
        if record.boiling_point.is_empty() {
            out_file.write_all(b"0.").unwrap();
            continue;
        }
        out_file.write_all(record.boiling_point.as_bytes()).unwrap();
        if record.boiling_point.find('.').is_none() {
            out_file.write_all(b".").unwrap();
        }
    }
    out_file.write_all(b"];\n").unwrap();
    out_file
        .write_all(b"const DENSITIES: [f32; 118] = [")
        .unwrap();
    first = true;
    for record in data.iter() {
        if first {
            first = false;
        } else {
            out_file.write_all(b", ").unwrap();
        }
        if record.density.is_empty() {
            out_file.write_all(b"0.").unwrap();
            continue;
        }
        out_file.write_all(record.density.as_bytes()).unwrap();
        if record.density.find('.').is_none() {
            out_file.write_all(b".").unwrap();
        }
    }
    out_file.write_all(b"];\n").unwrap();
    out_file
        .write_all(b"const GROUPS: [GroupBlock; 118] = [")
        .unwrap();
    first = true;
    for record in data.iter() {
        if first {
            first = false;
        } else {
            out_file.write_all(b", ").unwrap();
        }
        out_file
            .write_all(match &record.group_block.to_lowercase()[..] {
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
                p => panic!("Unknown group block: {}", p),
            })
            .unwrap();
    }
    out_file.write_all(b"];\n").unwrap();
    out_file
        .write_all(b"const YEARS_OF_DISCOVERED: [u16; 118] = [")
        .unwrap();
    first = true;
    for record in data.iter() {
        if first {
            first = false;
        } else {
            out_file.write_all(b", ").unwrap();
        }
        if record.year_discovered.is_empty() || record.year_discovered == "Ancient" {
            out_file.write_all(b"0").unwrap();
            continue;
        }
        out_file
            .write_all(record.year_discovered.as_bytes())
            .unwrap();
    }
    out_file.write_all(b"];\n").unwrap();
}
