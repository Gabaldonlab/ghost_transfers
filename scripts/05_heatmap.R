#!/usr/bin/env Rscript
#  Copyright Saioa Manzano-Morales, Moisès Bernabeu & Toni Gabaldón <saioa.manzano@bsc.es>
#  
#  This file is free software: you may copy, redistribute and/or modify it  
#  under the terms of the GNU General Public License as published by the  
#  Free Software Foundation, either version 2 of the License, or (at your  
#  option) any later version.  
#  
#  This file is distributed in the hope that it will be useful, but  
#  WITHOUT ANY WARRANTY; without even the implied warranty of  
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU  
#  General Public License for more details.  
#  
#  You should have received a copy of the GNU General Public License  
#  along with this program.  If not, see <http://www.gnu.org/licenses/>.  

results <- read.table("../outputs/donors_shift.tsv", sep="\t")

lookup <- data.frame(unique(c(results$donor1, results$donor2)))
colnames(lookup) <- "old" # Generate table with all donors
lookup$new <- ifelse(grepl(";", lookup$old), "Ancestor", lookup$old) # If composite, name "Ancestor"
# Put a letter next to each ancestor
indeces <- which(lookup$new == "Ancestor")
for (i in 1:length(indeces)) {
  lookup[indeces[i], "new"] <- paste0("Ancestor", LETTERS[i] )
}

results[,c(1,2)] <- lapply(results, function(x) lookup$new[match(x, lookup$old)]) # Change table according to lookup table
lookup$old <- gsub(";", "\n", lookup$old) # Aesthetic conversion
lookup <- lookup[grepl("Ancestor", lookup$new),c(2,1)] # Select ancestors
colnames(lookup) <- c("Ancestor", "Daughter phyla") #Change colnames

write.table(lookup, "../outputs/lookup.tsv", sep="\t", row.names=F)

colnames(results) <- c("donor1", "donor2", "shift_prop", "shift", "older_th_feca", "iter") # Colnames of results
write.table(data.frame(results), "../outputs/groups_results.tsv", sep="\t") # Store results

library(dplyr)
library(tidyr)
library(pheatmap)
library(viridis)

print("Plotting...")
# Pair is name of donors ordered alphabetically
results$pair <- ifelse(results$donor1 < results$donor2, paste0(results$donor1,";",
                                                               results$donor2),
                       paste0(results$donor2, ";", results$donor1))
results <- data.frame(results %>%
                        group_by(pair) %>%
                        summarise(n_positive = sum(shift), n_negative = sum(!shift))) #  Group by pair and count

results$prop <- results$n_positive/(results$n_negative + results$n_positive) # Calculate
results <- results %>% separate(pair, c('donor1', 'donor2')) # Separate pair names again

mat <- spread(results[,c(1,2,5)], donor2, prop) # Generate matrix
rownames(mat) <- mat[,1]; mat <- mat[,-1]; diag(mat) <- NA # Cleanup matrix

for (i in 1:dim(mat)[1]) {
  for (j in 1:dim(mat)[2]) {
    if (is.na(mat[i, j])) {
      mat[i, j] <- mat[j, i]
    }
  }
}

pheatmap(mat, cutree_rows = 4, cutree_cols = 4, na_col = 'white',
         filename = '../outputs/heatmap.pdf', width = 8, height = 8)

## Calc proportion of shifts older than FECA for pair Actinobacteria-Alphaproteobacteria ##
Alpha_actino <- results %>%
  group_by(pair, shift, older_th_feca) %>% filter(pair == "Actinobacteria;Proteobacteria") %>%  summarise(n = n()) %>%
  mutate(freq = n / sum(n))
                           
