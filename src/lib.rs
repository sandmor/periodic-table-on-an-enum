//! This is a complete periodic table for rust, including the following fields:
//! * Atomic number
//! * Symbol
//! * Name
//! * Atomic mass
//! * CPK Color
//! * Electron configuration
//! * Electronegativity
//! * AtomicRadius
//! * Ionization energy
//! * Electron affinity
//! * Oxidation states
//! * Standard state
//! * Melting point
//! * Boiling point
//! * Density
//! * Group block
//! * Year discovered

use std::mem;

include!(concat!(env!("OUT_DIR"), "/data.rs"));

#[derive(Debug, Copy, Clone, Eq, PartialEq, Ord, PartialOrd, Hash)]
pub enum GroupBlock {
    AlkaliMetal,
    AlkalineEarthMetal,
    Lanthanide,
    Actinide,
    TransitionMetal,
    PostTransitionMetal,
    Metalloid,
    NonMetal,
    Halogen,
    NobleGas,
}

macro_rules! lookup {
    ($table:expr, $term:expr) => {{
        let mut l = 0;
        let mut r = 117;
        while l <= r {
            let m = l + (r - l) / 2;
            if $table[m].0 == $term {
                return Some(unsafe { Element::from_id($table[m].1) });
            }
            if $table[m].0 < $term {
                l = m + 1;
            }
            if $table[m].0 > $term {
                r = m - 1;
            }
        }
        None
    }};
}

impl Element {
    pub fn get_oxidation_states(&self) -> &'static [i8] {
        if *self as u8 == 117 {
            // Last element of the list
            return &OXIDATION_STATES_DATA[OXIDATION_STATES[*self as usize].0 as usize..][..];
        }
        &OXIDATION_STATES_DATA[OXIDATION_STATES[*self as usize].0 as usize..
            OXIDATION_STATES[*self as usize].0 as usize+OXIDATION_STATES[*self as usize].1 as usize][..]
    }

    /// The id is the atomic number starting at zero
    #[inline(always)]
    pub unsafe fn from_id(id: u8) -> Element {
        mem::transmute(id)
    }

    pub fn from_symbol(sym: &str) -> Option<Element> {
        lookup!(SYMBOLS_SORTED_ALPHABETICALLY, sym)
    }


    /// Name must be lowercase
    pub fn from_name(name: &str) -> Option<Element> {
        lookup!(LOWERCASE_NAMES_SORTED_ALPHABETICALLY, name)
    }

    pub fn from_atomic_number(z: usize) -> Option<Element> {
        if z > 118 || z == 0 {
            return None;
        }
        Some(unsafe { mem::transmute((z-1) as u8) })
    }

    #[inline(always)]
    pub fn get_atomic_number(&self) -> usize {
        *self as usize + 1
    }

    #[inline(always)]
    pub fn get_atomic_mass(&self) -> f32 {
        ATOMIC_MASSES[*self as usize]
    }

    #[inline(always)]
    pub fn get_atomic_radius(&self) -> u16 {
        ATOMIC_RADIUS[*self as usize]
    }

    #[inline(always)]
    pub fn get_electronegativity(&self) -> f32 {
        ELECTRONEGATIVITIES[*self as usize]
    }

    #[inline(always)]
    pub fn get_electron_affinity(&self) -> f32 {
        ELECTRON_AFFINITIES[*self as usize]
    }

    #[inline(always)]
    pub fn get_electron_configuration(&self) -> &'static str {
        ELECTRON_CONFIGURATIONS[*self as usize]
    }

    #[inline(always)]
    pub fn get_ionization_energy(&self) -> f32 {
        IONIZATION_ENERGIES[*self as usize]
    }

    #[inline(always)]
    pub fn get_density(&self) -> f32 {
        DENSITIES[*self as usize]
    }

    #[inline(always)]
    pub fn get_melting_point(&self) -> f32 {
        MELTING_POINTS[*self as usize]
    }

    #[inline(always)]
    pub fn get_boiling_point(&self) -> f32 {
        BOILING_POINTS[*self as usize]
    }

    #[inline(always)]
    pub fn get_standard_state(&self) -> StateOfMatter {
        STANDARD_STATES[*self as usize]
    }

    #[inline(always)]
    pub fn get_symbol(&self) -> &'static str {
        SYMBOLS[*self as usize]
    }

    #[inline(always)]
    pub fn get_name(&self) -> &'static str {
        NAMES[*self as usize]
    }

    #[inline(always)]
    pub fn get_year_discovered(&self) -> u16 {
        YEARS_OF_DISCOVERED[*self as usize]
    }

    #[inline(always)]
    pub fn get_group(&self) -> GroupBlock {
        GROUPS[*self as usize]
    }

    #[inline(always)]
    pub fn get_cpk(&self) -> [u8; 3] {
        CPK[*self as usize]
    }

    /// The id is the atomic number starting at zero
    #[inline(always)]
    pub fn get_id(&self) -> u8 {
        *self as u8
    }
}
#[derive(Debug, Copy, Clone, Eq, PartialEq, Ord, PartialOrd, Hash)]
pub enum StateOfMatter {
    Solid,
    Liquid,
    Gas
}
