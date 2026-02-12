use extendr_api::prelude::*;
use gb_io::reader::SeqReader;
use gb_io::writer::SeqWriter;
use gb_io::seq::{self, Seq, Source, Topology}; 
use std::fs::File;
use std::borrow::Cow;

#[extendr]
fn write_genbank_rs(path: &str, sequence: &str, features: Dataframe<Robj>, metadata: List) -> Robj {
    // 1. Create Seq object
    let mut seq = Seq::empty();
    
    // Set sequence
    seq.seq = sequence.as_bytes().to_vec();
    
    // Set Metadata
    let get_str = |key: &str| -> Option<String> {
        metadata.dollar(key).ok().and_then(|x| x.as_str().map(|s| s.to_string()))
    };

    seq.definition = get_str("definition");
    seq.accession = get_str("accession");
    seq.version = get_str("version");
    
    // Name (LOCUS)
    seq.name = get_str("name");

    // Date
    if let Some(d_str) = get_str("date") {
        seq.date = parse_date(&d_str);
    }
    
    // Keywords
    if let Ok(k_robj) = metadata.dollar("keywords") {
        if let Some(strs) = k_robj.as_string_vector() {
             seq.keywords = Some(strs.join("; "));
        }
    }
    
    // Source
    if let Some(src_str) = get_str("source") {
        let source = Source {
             source: src_str,
             organism: get_str("organism"),
        };
        seq.source = Some(source);
    } 
    
    // Molecule Type (String in gb-io 0.9.0)
    if let Some(mol) = get_str("molecule_type") {
        seq.molecule_type = Some(mol);
    }
    
    // Division
    if let Some(div) = get_str("division") {
        seq.division = div;
    }
    
    // Topology
    if let Some(topo) = get_str("topology") {
         let topo_lower = topo.to_lowercase();
         if topo_lower == "linear" {
             seq.topology = Topology::Linear;
         } else if topo_lower == "circular" {
             seq.topology = Topology::Circular;
         }
    }
    
    // Features
    let keys_robj = features.dollar("key");
    let locs_robj = features.dollar("location");
    let quals_robj = features.dollar("qualifiers");
    
    if let (Ok(starts), Ok(locs), Ok(quals_list)) = (keys_robj, locs_robj, quals_robj) {
         let keys_vec: Vec<String> = starts.as_string_vector().unwrap_or_default();
         let locs_vec: Vec<String> = locs.as_string_vector().unwrap_or_default();
         
         if let Some(quals) = quals_list.as_list() {
             let n = keys_vec.len();
             if locs_vec.len() == n && quals.len() == n {
             
                 for i in 0..n {
                     let key_str = &keys_vec[i];
                     let loc_str = &locs_vec[i];
                     let qual_robj = quals.elt(i).unwrap(); 
                     
                     // FeatureKind is Cow<str>
                     let kind = Cow::from(key_str.clone());
                     
                     let location = match seq::Location::from_gb_format(loc_str) {
                         Ok(l) => l,
                         Err(e) => return Robj::from(format!("Error parsing location '{}': {}", loc_str, e))
                     };
                     
                     let mut f_quals: Vec<(Cow<str>, Option<String>)> = Vec::new(); 
                     
                     if let Some(q_strs) = qual_robj.as_string_vector() {
                          let names: Vec<String> = qual_robj.names().unwrap_or_default().map(|s| s.to_string()).collect(); 
                          
                          for (k, v) in names.iter().zip(q_strs.iter()) {
                              f_quals.push((Cow::from(k.to_string()), Some(v.clone())));
                          }
                     }
                     
                     let feature = gb_io::seq::Feature {
                         kind,
                         location,
                         qualifiers: f_quals,
                     };
                     
                     seq.features.push(feature);
                 }
             }
         }
    }

    // References
    if let Ok(refs_robj) = metadata.dollar("references") {
        if let Some(refs_list) = refs_robj.as_list() {
            for r_robj in refs_list.values() {
                if let Some(r_list) = r_robj.as_list() {
                    let get_ref_str = |k: &str| -> Option<String> {
                        r_list.dollar(k).ok().and_then(|x| x.as_str().map(|s| s.to_string()))
                    };
                    
                    let reference = seq::Reference {
                        description: get_ref_str("description").unwrap_or_default(),
                        authors: get_ref_str("authors"),
                        consortium: get_ref_str("consortium"),
                        title: get_ref_str("title").unwrap_or_default(),
                        journal: get_ref_str("journal"),
                        pubmed: get_ref_str("pubmed"),
                        remark: get_ref_str("remark"),
                    };
                    seq.references.push(reference);
                }
            }
        }
    }

    // Write file
    let file = match File::create(path) {
        Ok(f) => f,
        Err(e) => return Robj::from(format!("Error creating file: {}", e))
    };
    let mut writer = SeqWriter::new(file);
    if let Err(e) = writer.write(&seq) {
        return Robj::from(format!("Error writing genbank: {}", e));
    }

    Robj::from(true)
}

fn parse_date(s: &str) -> Option<seq::Date> {
    let parts: Vec<&str> = s.split('-').collect();
    if parts.len() != 3 { return None; }
    
    let day = parts[0].parse::<u32>().ok()?;
    let month_str = parts[1].to_uppercase();
    let month = match month_str.as_str() {
        "JAN" => 1, "FEB" => 2, "MAR" => 3, "APR" => 4,
        "MAY" => 5, "JUN" => 6, "JUL" => 7, "AUG" => 8,
        "SEP" => 9, "OCT" => 10, "NOV" => 11, "DEC" => 12,
        _ => return None
    };
    let year = parts[2].parse::<i32>().ok()?;
    
    seq::Date::from_ymd(year, month, day).ok()
}

#[extendr]
fn read_genbank_rs(path: &str) -> Robj {
    let file_result = File::open(path);
    if let Err(e) = file_result {
        return Robj::from(format!("Error opening file: {}", e));
    }
    let file = file_result.unwrap();
    let reader = SeqReader::new(file);
    
    let mut records = Vec::new();
    
    for seq_res in reader {
        match seq_res {
            Ok(seq) => {
                // Metadata
                let keywords_vec = match &seq.keywords {
                     Some(k) => k.split("; ").map(|s| s.to_string()).collect(),
                     None => Vec::new()
                };
                
                // References
                let mut refs_list = Vec::new();
                for r in seq.references {
                    let r_list = list!(
                        description = r.description,
                        authors = r.authors,
                        consortium = r.consortium,
                        title = r.title,
                        journal = r.journal,
                        pubmed = r.pubmed,
                        remark = r.remark
                    );
                    refs_list.push(r_list);
                }
                
                let metadata = list!(
                    name = seq.name.clone().unwrap_or_default(),
                    definition = seq.definition.clone().unwrap_or_default(),
                    accession = seq.accession.clone().unwrap_or_default(),
                    version = seq.version.clone().unwrap_or_default(),
                    keywords = keywords_vec,
                    source = seq.source.as_ref().map(|s| s.source.clone()).unwrap_or_default(),
                    organism = seq.source.as_ref().map(|s| s.organism.clone()).unwrap_or_default(),
                    molecule_type = seq.molecule_type.clone().unwrap_or_default(),
                    topology = seq.topology.to_string(),
                    division = seq.division.to_string(),
                    date = seq.date.as_ref().map(|d| d.to_string()),
                    references = List::from_values(refs_list)
                );

                // Features
                let mut kinds = Vec::with_capacity(seq.features.len());
                let mut locs = Vec::with_capacity(seq.features.len());
                let mut quals_col = Vec::with_capacity(seq.features.len());
                
                for f in seq.features {
                    kinds.push(f.kind.to_string());
                    locs.push(f.location.to_string());
                    
                    let mut q_names = Vec::with_capacity(f.qualifiers.len());
                    let mut q_vals = Vec::with_capacity(f.qualifiers.len());
                    
                    for (k, v) in f.qualifiers {
                         q_names.push(k.to_string()); 
                         q_vals.push(v.unwrap_or_default());
                    }
                     
                    let mut q_robj = Robj::from(q_vals);
                    if !q_names.is_empty() {
                         let _ = q_robj.set_names(&q_names);
                    }
                    quals_col.push(q_robj);
                }
                
                let features_list = list!(
                    key = kinds,
                    location = locs,
                    qualifiers = List::from_values(quals_col)
                );
                
                // Sequence
                let sequence = String::from_utf8_lossy(&seq.seq).to_string();
                
                let record = list!(
                    metadata = metadata,
                    features = features_list,
                    sequence = sequence
                );
                
                records.push(record);
            },
            Err(e) => return Robj::from(format!("Error reading sequence: {}", e))
        }
    }
    
    // Return list of records
    List::from_values(records).into()
}

extendr_module! {
    mod rgbio;
    fn read_genbank_rs;
    fn write_genbank_rs;
}
