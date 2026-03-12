library(ape)
library(dplyr)

# 1. Carica l'albero
tree <- read.tree("tree.nwk")

# 2. Carica la tabella di tip e ordini
tax <- read.table("taxonomy.tsv", header = TRUE, sep="\t", stringsAsFactors = FALSE)

# 3. Prepara un data frame vuoto per i risultati
results <- data.frame(order = character(),
                      is_monophyletic = logical(),
                      bootstrap = numeric(),
                      stringsAsFactors = FALSE)

# 4. Controllo per ogni ordine
for(ord in unique(tax$order)){
  
  tips_in_order <- tax$tip[tax$order == ord]
  
  # Controlla se sono monofiletici
  mono <- is.monophyletic(tree, tips_in_order)
  
  # Se monofiletico, trova il nodo MRCA
  if(mono){
    mrca_node <- getMRCA(tree, tips_in_order)
    
    # Bootstrap del nodo (se disponibile)
    if(!is.null(tree$node.label) && length(tree$node.label) >= mrca_node - length(tree$tip.label)){
      # Nodo interno i numeri da Ntip+1 a Ntip+Nnode
      bootstrap_value <- as.numeric(tree$node.label[mrca_node - length(tree$tip.label)])
    } else {
      bootstrap_value <- NA
    }
    
  } else {
    bootstrap_value <- NA
  }
  
  # Aggiungi al data frame
  results <- rbind(results, data.frame(order = ord,
                                       is_monophyletic = mono,
                                       bootstrap = bootstrap_value,
                                       stringsAsFactors = FALSE))
}

# 5. Mostra i risultati
print(results)