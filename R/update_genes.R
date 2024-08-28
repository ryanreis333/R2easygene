#' Update Gene Symbols and Normalize Seurat Object
#'
#' This function updates gene names in a Seurat object using the latest gene symbol mapping
#' table from the HGNChelper package. It filters out low-expressed genes, updates gene names
#' to their approved symbols, and handles duplicated genes. Optionally, it can return
#' a new Seurat object with updated metadata, including mitochondrial and ribosomal gene
#' percentages.
#'
#' @param Seurat_obj A Seurat object. The RNA or any specified assay counts will be used.
#' @param assay The name of the assay from which to extract the counts matrix. Default is "RNA".
#' @param return_seurat Logical. If TRUE, a new Seurat object is returned with the updated gene names
#' and percentage features for mitochondrial and ribosomal genes. If FALSE, the function returns the
#' updated counts matrix. Default is TRUE.
#' @param min.cells Numeric. The minimum number of cells in which a gene must be detected to be
#' retained in the Seurat object. Default is 3.
#' @param min.features Numeric. The minimum number of features (genes) that must be detected in a
#' cell for that cell to be retained in the Seurat object. Default is 200.
#'
#' @return If \code{return_seurat} is TRUE, the function returns a new Seurat object
#' with updated gene names and two additional metadata columns: \code{percent.mt}, which indicates
#' the percentage of mitochondrial genes, and \code{percent.rb}, which shows the percentage
#' of ribosomal genes. If \code{return_seurat} is FALSE, the function returns the updated
#' counts matrix with corrected gene names.
#'
#' @details This function leverages the HGNChelper package to ensure that all gene names
#' in a Seurat object match the currently approved HGNC symbols. It filters out genes
#' detected in fewer than \code{min.cells} cells and performs several operations to
#' handle duplicated gene symbols. For downstream analysis, it adds the proportion of
#' mitochondrial genes (if applicable) and ribosomal genes, accounting for both human
#' and mouse species.
#'
#' @examples
#' \dontrun{
#' # Example with a Seurat object
#' updated_seurat <- update_genes(Seurat_obj, assay = "RNA")
#' updated_matrix <- update_genes(Seurat_obj, assay = "RNA", return_seurat = FALSE)
#' }
#'
#' @import Seurat
#' @importFrom HGNChelper hgnc.table
#' @importFrom Seurat PercentageFeatureSet CreateSeuratObject
#' @seealso \code{\link[Seurat]{CreateSeuratObject}}, \code{\link[Seurat]{PercentageFeatureSet}}, \code{\link[HGNChelper]{hgnc.table}}
#' @export
#'

update_genes <- function(Seurat_obj, assay = "RNA", return_seurat = TRUE, min.cells = 3, min.features = 200) {

# Load gene symbol mapping table from HGNChelper
hgnc.table <- HGNChelper::hgnc.table

# Load gene symbols
all_genes <- unique(c(hgnc.table$Symbol, hgnc.table$Approved.Symbol))
approved_genes <- unique(hgnc.table$Approved.Symbol)

# Extract counts matrix from the specified assay in Seurat object
mat <- Seurat_obj@assays$assay@counts

# Check if the counts matrix is valid
if (!is.matrix(mat) && !is.data.frame(mat)) {
  stop("The counts matrix is not valid or does not exist.")
}

# Filter genes with low counts and non-matching genes
keep_genes <- rownames(mat)[rowSums(mat) >= min.cells & rownames(mat) %in% all_genes]
mat <- mat[keep_genes, ]

# Fix gene symbols
symbol_map <- setNames(hgnc.table$Approved.Symbol, hgnc.table$Symbol)
new_names <- rownames(mat)
mapped_indices <- match(new_names, names(symbol_map))
valid_indices <- !is.na(mapped_indices)
new_names[valid_indices] <- symbol_map[mapped_indices[valid_indices]]

# Sort both mat and new_names
order_idx <- order(new_names)
mat <- mat[order_idx, ]
new_names <- new_names[order_idx]

# Identify duplicate genes
duplicate_genes <- new_names[duplicated(new_names) | duplicated(new_names, fromLast = TRUE)]
unique_genes <- setdiff(new_names, duplicate_genes)

# Subset the matrix
unique_mat <- mat[new_names %in% unique_genes, , drop = FALSE]
duplicate_mat <- mat[new_names %in% duplicate_genes, , drop = FALSE]

# Sum the duplicate rows
summed_duplicates <- rowsum(duplicate_mat, new_names[new_names %in% duplicate_genes])

# Combine the results
final_mat <- rbind(unique_mat, summed_duplicates)

# Sort the final matrix by gene names
final_mat <- final_mat[order(rownames(final_mat)), ]
rownames(final_mat) <- symbol_map[rownames(final_mat)]

if (return_seurat) {
  # Create a new Seurat object using the updated matrix
  Seurat_obj <- CreateSeuratObject(counts = final_mat, meta.data = Seurat_obj[[]], min.cells = min.cells, min.features = min.features)

  # Add percentage of mitochondrial genes
  if (any(grepl("^MT-", rownames(final_mat), ignore.case = TRUE))) {
    Seurat_obj[["percent.mt"]] <- PercentageFeatureSet(Seurat_obj, pattern = "^MT-", assay = assay)
  } else {
    Seurat_obj[["percent.mt"]] <- PercentageFeatureSet(Seurat_obj, pattern = "^mt-", assay = assay)
  }

  # Add percentage of ribosomal genes
  Seurat_obj[["percent.rb"]] <- PercentageFeatureSet(Seurat_obj, pattern = "^RPS|^RPL|^Rpl|^Rps", assay = assay)

  return(Seurat_obj)
} else {
  return(final_mat)
}
}

