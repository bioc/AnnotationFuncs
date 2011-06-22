#' Performs quicker lookup for orthologs in homologe data packages
#'
#' Using the INPARANOID data packages such as \code{hom.Hs.inp.db} is very, very slow and can take up to 11 min (on this particular developers workstation).
#' This function introduces a new method that can do it in just 20 seconds (on the developers workstation).
#' In addition, it includes options for translating between different identifers both before and after the mapping.
#' 
#' @param values Vector, coerced to character vector, of values needed mapping by homology.
#' @param mapping Homology mapping object, such as \code{hom.Hs.inpBOSTA} or \code{revmap(hom.Hs.inpBOSTA)}.
#' @param genus Character vector. 5 character INPARANOID style genus name of the mapping object, e.g. 'BOSTA' for both \code{hom.Hs.inpBOSTA} and \code{revmap(hom.Hs.inpBOSTA)}.
#' @param threshold Numeric value between 0 and 1. Only clustered homologues with a parwise score above the threshold is included.
#'                  The native implementation has this set to 1.
#' @param pre.from Mapping object if \code{values} needs translation before mapping. 
#'             E.g. \code{values} are entrez and \code{hom.Hs.inpBOSTA} requires ENSEMBLPROT, \code{hom.Hs.inpAPIME} requires Refseq (?).
#'             Arguments \code{from} and \code{to} are just like in \code{\link{translate}}.
#' @param pre.to Second part of translation before mapping.
#' @param post.from Translate the result from homology mapping to a desired id; just like in \code{\link{translate}}.
#' @param post.to Second part of translation after mapping.
#' @param ... Additional arguments sent to \code{\link{translate}}.
#' @return List. Names of list corresponds to \code{values}, except those that could not be mapped nor translated.
#'               Entries are character vectors.
#' @references \code{?hom.Hs.inp.db} - \url{http://inparanoid.sbc.su.se/}
#' 
#'  Berglund, A.C., Sjolund, E., Ostlund, G., Sonnhammer, E.L.L. (2008) 
#'  InParanoid 6: eukaryotic ortholog clusters with inparalogs
#'  \emph{Nucleic Acids Res.} \bold{36}:D263--266
#' 
#'  O'Brien, K.P., Maido, R., Sonnhammer, E.L.L (2005)
#'  Inparanoid: A Comprehensive Database of Eukaryotic Orthologs
#'  \emph{NAR} \bold{33}:D476--D480
#'
#'  Remm, M., Storm, C.E.V, Sonnhammer, E.L.L (2001)
#'  Automatic clustering of orthologs and in-paralogs from pairwise species comparisons
#'  \emph{J. Mol. Biol.} \bold{314}:1041--1052
#'  
#' @seealso \code{\link{translate}}, \code{\link{.getTableName}}, \code{\link{mapLists}}
#' @export
#' @author Stefan McKinnon Edwards \email{stefan.hoj-edwards@@agrsci.dk}
#' @examples
#' library(hom.Hs.inp.db)
#' library(org.Hs.eg.db)
#' library(org.Bt.eg.db)
#' getOrthologs("ENSBTAP00000024572", revmap(hom.Hs.inpBOSTA), 'BOSTA') 
#' # And now, we will map from entrez genes 1, 2 and 3 to bovine Refseq
#' bovine.ensembl <- getOrthologs(c(1,2,3), hom.Hs.inpBOSTA, 'BOSTA', pre.from=org.Hs.egENSEMBLPROT, post.from=org.Bt.egENSEMBLPROT2EG)
#' refseqs <- translate(unlist(bovine.ensembl, use.names=FALSE), org.Bt.egREFSEQ)
#' hs2bt.refseqs <- mapLists(bovine.ensembl, refseqs)
#' # Another way of doing it:
#' hs2bt.refseqs2 <- lapply(bovine.ensembl, translate, from=org.Bt.egREFSEQ, simplify=TRUE) # simplify=TRUE is very important here!

#library(AnnotationFuncs)
#library(org.Bt.eg.db)
#library(hom.Hs.inp.db)
#setwd('C:/TXT/Joanna')
#load('genelist.Rdata')
#ens <- translate(genes, org.Bt.egENSEMBLPROT)
#values <- unlist(ens, use.names=F)
#mapping <- hom.Hs.inpBOSTA
#genus <- 'BOSTA'
#threshold <- 1
#tbl <- 'Bos_taurus'

getOrthologs <- function(values, mapping, genus, threshold=1, 
		pre.from=NULL, pre.to=NULL, 
		post.from=NULL, post.to=NULL, 
		...) {
    values <- as.character(values)
    values <- unique(values)
    
    threshold <- as.numeric(threshold)
    threshold <- max(0, min(threshold, 1))
        
    # Check if we do some translating first:
    if (!is.null(pre.from)) {
        trans1 <- translate(values, pre.from, pre.to, ...)
        values <- unlist(trans1, use.names=F)
    }
    
    # Check validity of `mapping` and genus.
    stopifnot(dbmeta(dbconn(mapping), "DBSCHEMA") == 'INPARANOID_DB')
    stopifnot(nchar(genus) == 5)
    .dbEscapeString(genus)
    genus <- toupper(genus)
    
    # Get the table name and check that it exists
    tbl <- .getTableName(genus)
    
    if ((length(tbl) == 0) | !(tbl %in% dbListTables(dbconn(mapping))))      
        stop('Provided genus was not recognized for the mapping object.')
    
    # Do the DB magic
    conn <- dbConnect(dbDriver('SQLite'), '')   # Create a temporary table.
    dbSendQuery(conn, paste("ATTACH DATABASE '",dbfile(mapping),"' AS hom", sep=''))
    res <- dbSendQuery(conn, 'CREATE TABLE input (id1 TEXT);')
    dbClearResult(res)
    dbWriteTable(conn, 'input', data.frame(id1=values), row.names = FALSE, overwrite=FALSE, append=TRUE )    
    #dbSendQuery(conn, 'BEGIN TRANSACTION')
    #bwahh <- sapply(values, function(s) dbSendQuery(conn, paste('INSERT INTO input (id1) VALUES ("',s,'")', sep='')))
    #dbSendQuery(conn, 'COMMIT TRANSACTION')

    # Proceed as normal, now that we have loaded our data. :)
    dbSendQuery(conn, sprintf('CREATE TABLE clusters AS SELECT id1, clust_id FROM input INNER JOIN hom.%s as inp ON inp.inp_id=id1 WHERE inp.score >= %i AND inp.seed_status = "100%%"', tbl, threshold))
    #res <- dbGetQuery(conn, sprintf('SELECT id1, clust_id FROM input INNER JOIN hom.%s as inp ON inp.inp_id=id1', tbl))
    isnot <- ifelse(direction(mapping) == 1, 'IS', 'IS NOT')
    res <- dbGetQuery(conn, sprintf('SELECT id1, inp_id FROM clusters INNER JOIN hom.%s as inp ON inp.clust_id=clusters.clust_id WHERE inp.score >= %i AND inp.species %s "%s" AND inp.seed_status = "100%%"', tbl, threshold, isnot, genus))
    
    dbDisconnect(conn)
    
    # A bit of cleaning up
    #res.list <- unstack(res, res$inp_id ~ res$id1)
    res.list <- split(res$inp_id, f=factor(res$id1, levels=unique(res$id1)))
    
    if (!is.null(pre.from)) {
        res.list <- mapLists(trans1, res.list)
    }
    if (!is.null(post.from)) {
        trans2 <- translate(res$inp_id, post.from, post.to, ...)
        res.list <- mapLists(res.list, trans2)
    }
    res.list <- lapply(res.list, unique)
    
    return(res.list)
}


#'  Gets the table name from the INPARANOID style genus names.
#'
#'  The INPARANOID style genus name is a 5 letter acronym of the species name. 
#'  Quote INPARANOID (\code{?hom.Hs.inpBOSTA}):
#'
#'  \emph{Names for these maps are done in the "INPARANOID style" which means that they are normally the 1st three letters of the genus followed by the 1st two letters of the species. For example: "Mus musculus" becomes "MUSMU", "Homo sapiens" becomes "HOMSA", "Monodelphis domestica" becomes "MONDO" etc. This means that for most of these organisms it will be possible to easily guess the abbreviations used. An exception may occur in the future if a new model organism has a very similar genus and species name to an existing one. }
#'
#' 
#'  @param genus 5 character INPARANOID genus name, such as "BOSTA", "HOMSA" or "MUSMU".
#'  @return Table name for genus.
#'  @references \url{http://www.bioconductor.org/packages/release/bioc/html/AnnotationDbi.html}
#'  @author Stefan McKinnon Edwards \email{stefanm.edwards@@agrsci.dk}
.getTableName <- function(genus) {
    # Find the AnnotationDbi source
    # Open createAnnObjs.INPARANOID_DB.R
    # Find the long bit of code that create variable 'fields' with all the species and the 5-character genus names.
    # Create variable fields from it and use the next three lines to swap names and values.
    # tmp <- fields
    # fields <- names(fields)
    # names(fields) <- tmp
    # Execute 'fix(fields)' and paste the result down here:
    fields <- structure(c("Acyrthosiphon_pisum", "Aedes_aegypti", "Anopheles_gambiae", 
    "Apis_mellifera", "Arabidopsis_thaliana", "Aspergillus_fumigatus", 
    "Batrachochytrium_dendrobatidis", "Bombyx_mori", "Bos_taurus", 
    "Branchiostoma_floridae", "Brugia_malayi", "Caenorhabditis_brenneri", 
    "Caenorhabditis_briggsae", "Caenorhabditis_elegans", "Caenorhabditis_japonica", 
    "Caenorhabditis_remanei", "Candida_albicans", "Candida_glabrata", 
    "Canis_familiaris", "Capitella_spI", "Cavia_porcellus", "Chlamydomonas_reinhardtii", 
    "Ciona_intestinalis", "Ciona_savignyi", "Coccidioides_immitis", 
    "Coprinopsis_cinereus", "Cryptococcus_neoformans", "Cryptosporidium_hominis", 
    "Cryptosporidium_parvum", "Culex_pipiens", "Cyanidioschyzon_merolae", 
    "Danio_rerio", "Daphnia_pulex", "Debaryomyces_hansenii", "Dictyostelium_discoideum", 
    "Drosophila_ananassae", "Drosophila_grimshawi", "Drosophila_melanogaster", 
    "Drosophila_mojavensis", "Drosophila_pseudoobscura", "Drosophila_virilis", 
    "Drosophila_willistoni", "Entamoeba_histolytica", "Equus_caballus", 
    "Escherichia_coliK12", "Fusarium_graminearum", "Gallus_gallus", 
    "Gasterosteus_aculeatus", "Giardia_lamblia", "Helobdella_robusta", 
    "Homo_sapiens", "Ixodes_scapularis", "Kluyveromyces_lactis", 
    "Leishmania_major", "Lottia_gigantea", "Macaca_mulatta", "Magnaporthe_grisea", 
    "Monodelphis_domestica", "Monosiga_brevicollis", "Mus_musculus", 
    "Nasonia_vitripennis", "Nematostella_vectensis", "Neurospora_crassa", 
    "Ornithorhynchus_anatinus", "Oryza_sativa", "Oryzias_latipes", 
    "Ostreococcus_tauri", "Pan_troglodytes", "Pediculus_humanus", 
    "Physcomitrella_patens", "Phytophthora_ramorum", "Phytophthora_sojae", 
    "Plasmodium_falciparum", "Plasmodium_vivax", "Pongo_pygmaeus", 
    "Populus_trichocarpa", "Pristionchus_pacificus", "Puccinia_graminis", 
    "Rattus_norvegicus", "Rhizopus_oryzae", "Saccharomyces_cerevisiae", 
    "Schistosoma_mansoni", "Schizosaccharomyces_pombe", "Sclerotinia_sclerotiorum", 
    "Sorghum_bicolor", "Stagonospora_nodorum", "Strongylocentrotus_purpuratus", 
    "Takifugu_rubripes", "Tetrahymena_thermophila", "Tetraodon_nigroviridis", 
    "Thalassiosira_pseudonana", "Theileria_annulata", "Theileria_parva", 
    "Tribolium_castaneum", "Trichomonas_vaginalis", "Trichoplax_adhaerens", 
    "Trypanosoma_cruzi", "Ustilago_maydis", "Xenopus_tropicalis", 
    "Yarrowia_lipolytica"), .Names = c("ACYPI", "AEDAE", "ANOGA", 
    "APIME", "ARATH", "ASPFU", "BATDE", "BOMMO", "BOSTA", "BRAFL", 
    "BRUMA", "CAEBRE", "CAEBR", "CAEEL", "CAEJA", "CAERE", "CANAL", 
    "CANGL", "CANFA", "CAPSP", "CAVPO", "CHLRE", "CIOIN", "CIOSA", 
    "COCIM", "COPCI", "CRYNE", "CRYHO", "CRYPA", "CULPI", "CYAME", 
    "DANRE", "DAPPU", "DEBHA", "DICDI", "DROAN", "DROGR", "DROME", 
    "DROMO", "DROPS", "DROVI", "DROWI", "ENTHI", "EQUCA", "ESCCO", 
    "FUSGR", "GALGA", "GASAC", "GIALA", "HELRO", "HOMSA", "IXOSC", 
    "KLULA", "LEIMA", "LOTGI", "MACMU", "MAGGR", "MONDO", "MONBR", 
    "MUSMU", "NASVI", "NEMVE", "NEUCR", "ORNAN", "ORYSA", "ORYLA", 
    "OSTTA", "PANTR", "PEDPA", "PHYPA", "PHYRA", "PHYSO", "PLAFA", 
    "PLAVI", "PONPY", "POPTR", "PRIPA", "PUCGR", "RATNO", "RHIOR", 
    "SACCE", "SCHMA", "SCHPO", "SCLSC", "SORBI", "STANO", "STRPU", 
    "TAKRU", "TETTH", "TETNI", "THAPS", "THEAN", "THEPA", "TRICA", 
    "TRIVA", "TRIAD", "TRYCR", "USTMA", "XENTR", "YARLI"))

    genus <- toupper(genus)
    return(removeNAs(fields[genus]))
}

#'  Private Escape string
#' 
#'  Does not escape strings, but raises an error if any character expect normal letters and underscores are found in the string.
#'  @param str String to test
#'  @param raise.error Logical, whether to raise an error or not.
#'  @return Invisible logical
#'  
.dbEscapeString <- function(str, raise.error=TRUE) {
    matches <- grepl('[^A-Za-z0-9_]+', str)
    res <- any(matches)
    if (raise.error & res) stop('Supplied SQL string contains illegal characters!')
    invisible(res)
}
