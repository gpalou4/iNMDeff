#library("Seurat")
library("httr")

# The following functions are extracted from the Seurat R package
# But I didn't want to install de whole package only for this

####
#' Get updated synonyms for gene symbols
#'
#' Find current gene symbols based on old or alias symbols using the gene
#' names database from the HUGO Gene Nomenclature Committee (HGNC)
#'
#' @details For each symbol passed, we query the HGNC gene names database for
#' current symbols that have the provided symbol as either an alias
#' (\code{alias_symbol}) or old (\code{prev_symbol}) symbol. All other queries
#' are \strong{not} supported.
#'
#' @note This function requires internet access
#'
#' @param symbols A vector of gene symbols
#' @param timeout Time to wait before canceling query in seconds
#' @param several.ok Allow several current gene symbols for each
#' provided symbol
#' @param search.types Type of query to perform:
#' \describe{
#'  \item{\dQuote{\code{alias_symbol}}}{Find alternate symbols for the genes
#'  described by \code{symbols}}
#'  \item{\dQuote{\code{prev_symbol}}}{Find new new symbols for the genes
#'  described by \code{symbols}}
#' }
#' This parameter accepts multiple options and short-hand options
#' (eg. \dQuote{\code{prev}} for \dQuote{\code{prev_symbol}})
#' @param verbose Show a progress bar depicting search progress
#' @param ... Extra parameters passed to \code{\link[httr]{GET}}
#'
#' @return \code{GeneSymbolThesarus}:, if \code{several.ok}, a named list
#' where each entry is the current symbol found for each symbol provided and
#' the names are the provided symbols. Otherwise, a named vector with the
#' same information.
#'
#' @source \url{https://www.genenames.org/} \url{https://www.genenames.org/help/rest/}
#'
#' @importFrom utils txtProgressBar setTxtProgressBar
#' @importFrom httr GET accept_json timeout status_code content
#'
#' @rdname UpdateSymbolList
#' @name UpdateSymbolList
#'
#' @export
#' @concept utilities
#'
#' @seealso \code{\link[httr]{GET}}
#'
#' @examples
#' \dontrun{
#' GeneSybmolThesarus(symbols = c("FAM64A"))
#' }
#'


GeneSymbolThesarus <- function(
  symbols,
  timeout = 10,
  several.ok = FALSE,
  search.types = c('alias_symbol', 'prev_symbol'),
  verbose = TRUE
) {
  db.url <- 'http://rest.genenames.org/fetch'
  # search.types <- c('alias_symbol', 'prev_symbol')
  search.types <- match.arg(arg = search.types, several.ok = TRUE)
  synonyms <- vector(mode = 'list', length = length(x = symbols))
  not.found <- vector(mode = 'logical', length = length(x = symbols))
  multiple.found <- vector(mode = 'logical', length = length(x = symbols))
  names(x = multiple.found) <- names(x = not.found) <- names(x = synonyms) <- symbols
  if (verbose) {
    pb <- txtProgressBar(max = length(x = symbols), style = 3, file = stderr())
  }
  for (symbol in symbols) {
    sym.syn <- character()
    for (type in search.types) {
      response <- GET(
        url = paste(db.url, type, symbol, sep = '/'),
        config = c(accept_json(), timeout(seconds = timeout))
      )
      if (!identical(x = status_code(x = response), y = 200L)) {
        next
      }
      response <- content(x = response)
      if (response$response$numFound != 1) {
        if (response$response$numFound > 1) {
          warning(
            "Multiple hits found for ",
            symbol,
            " as ",
            type,
            ", skipping",
            call. = FALSE,
            immediate. = TRUE
          )
        }
        next
      }
      sym.syn <- c(sym.syn, response$response$docs[[1]]$symbol)
    }
    not.found[symbol] <- length(x = sym.syn) < 1
    multiple.found[symbol] <- length(x = sym.syn) > 1
    if (length(x = sym.syn) == 1 || (length(x = sym.syn) > 1 && several.ok)) {
      synonyms[[symbol]] <- sym.syn
    }
    if (verbose) {
      setTxtProgressBar(pb = pb, value = pb$getVal() + 1)
    }
  }
  if (verbose) {
    close(con = pb)
  }
  if (sum(not.found) > 0) {
    warning(
      "The following symbols had no synonyms: ",
      paste(names(x = which(x = not.found)), collapse = ', '),
      call. = FALSE,
      immediate. = TRUE
    )
  }
  if (sum(multiple.found) > 0) {
    msg <- paste(
      "The following symbols had multiple synonyms:",
      paste(names(x = which(x = multiple.found)), sep = ', ')
    )
    if (several.ok) {
      message(msg)
      message("Including anyways")
    } else {
      warning(msg, call. = FALSE, immediate. = TRUE)
    }
  }
  synonyms <- Filter(f = Negate(f = is.null), x = synonyms)
  if (!several.ok) {
    synonyms <- unlist(x = synonyms)
  }
  return(synonyms)
}

#' @rdname UpdateSymbolList
#'
#' @return \code{UpdateSymbolList}: \code{symbols} with updated symbols from
#' HGNC's gene names database
#'
#' @export
#' @concept utilities
#'
#' @examples
#' \dontrun{
#' UpdateSymbolList(symbols = cc.genes$s.genes)
#' }
#'
UpdateSymbolList <- function(
  symbols,
  timeout = 10,
  several.ok = FALSE,
  verbose = TRUE,
  ...
) {
  new.symbols <- suppressWarnings(expr = GeneSymbolThesarus(
    symbols = symbols,
    timeout = timeout,
    several.ok = several.ok,
    search.types = 'prev_symbol',
    verbose = verbose,
    ...
  ))
  if (length(x = new.symbols) < 1) {
    warning("No updated symbols found", call. = FALSE, immediate. = TRUE)
  } else {
    if (verbose) {
      message("Found updated symbols for ", length(x = new.symbols), " symbols")
      x <- sapply(X = new.symbols, FUN = paste, collapse = ', ')
      message(paste(names(x = x), x, sep = ' -> ', collapse = '\n'))
    }
    for (sym in names(x = new.symbols)) {
      index <- which(x = symbols == sym)
      symbols <- append(
        x = symbols[-index],
        values = new.symbols[[sym]],
        after = index - 1
      )
    }
  }
  return(symbols)
}

####

# 1) Update Gene Symbols using UpdateSymbolList function from Seurat package

### TEMP ###
# Open ensembl GTF
library("rtracklayer")
gtf <- rtracklayer::import('/g/strcombio/fsupek_home/gpalou/data/NMD_project/conversor_tables/ENSEMBL/gencode.v26.annotation.gtf')
gtf_df <- as.data.frame(gtf)

# Symbols shared before updating
table(UPF1.dir.unconf.gene.symbols%in%unique(gtf_df$gene_name))
# Symbols shared after updating
table(UPF1.dir.unconf.gene.symbols.updated%in%unique(gtf_df$gene_name))

############

# 1.1) Tani (2012)

NMD_targets_path = "/g/strcombio/fsupek_home/gpalou/data/NMD_project/NMD_targets/"

UPF1.dir.conf.gene.symbols <- as.character(read.table(file = paste(NMD_targets_path, "raw/Tani_UPF1_dir_conf_gene_symbols.txt", sep = ""), 
                                         header = FALSE, sep = "\t")$V1)
UPF1.dir.conf.gene.symbols.updated <- unique(UpdateSymbolList(symbols = UPF1.dir.conf.gene.symbols))
write.table(UPF1.dir.conf.gene.symbols.updated , file = paste0(NMD_targets_path,"updated/Tani_UPF1_dir_conf_gene_symbols_updated.txt"), quote = FALSE, row.names = FALSE, col.names = FALSE)

UPF1.dir.unconf.gene.symbols <- as.character(read.table(file = paste(NMD_targets_path, "raw/Tani_UPF1_dir_unconf_gene_symbols.txt", sep = ""), 
                                         header = FALSE, sep = "\t")$V1)
UPF1.dir.unconf.gene.symbols.updated <- unique(UpdateSymbolList(symbols = UPF1.dir.unconf.gene.symbols))
write.table(UPF1.dir.unconf.gene.symbols.updated , file = paste0(NMD_targets_path, "updated/Tani_UPF1_dir_unconf_gene_symbols_updated.txt"), quote = FALSE, row.names = FALSE, col.names = FALSE)

UPF1.ind.conf.gene.symbols <- as.character(read.table(file = paste(NMD_targets_path, "raw/Tani_UPF1_ind_conf_gene_symbols.txt", sep = ""), 
                                         header = FALSE, sep = "\t")$V1)
UPF1.ind.conf.gene.symbols.updated <- unique(UpdateSymbolList(symbols = UPF1.ind.conf.gene.symbols))
write.table(UPF1.ind.conf.gene.symbols.updated , file = paste0(NMD_targets_path,"updated/Tani_UPF1_ind_conf_gene_symbols_updated.txt"), quote = FALSE, row.names = FALSE, col.names = FALSE)

# UPF1 all
UPF1.all.gene.symbols.updated <- unique(c(UPF1.dir.conf.gene.symbols.updated,UPF1.dir.unconf.gene.symbols.updated,UPF1.ind.conf.gene.symbols.updated))
write.table(UPF1.all.gene.symbols.updated , file = paste0(NMD_targets_path,"updated/Tani_UPF1_all_gene_symbols_updated.txt"), quote = FALSE, row.names = FALSE, col.names = FALSE)

# 1.2) Schmidt (2015)
SMG6.gene.symbols <- as.character(read.table(file = paste(NMD_targets_path, "raw/Schmidt_SMG6_gene_symbol.txt", sep = ""), header = FALSE, sep = "\t")$V1)
SMG6.gene.symbols.updated <- unique(UpdateSymbolList(symbols = SMG6.gene.symbols))
write.table(SMG6.gene.symbols.updated, file = paste0(NMD_targets_path,"updated/Schmidt_SMG6_gene_symbol_updated.txt"), quote = FALSE, row.names = FALSE, col.names = FALSE)

