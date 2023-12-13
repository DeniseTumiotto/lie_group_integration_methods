#include "expandconfig.h"

#ifdef DEBUG
#define DEBUG_MSG_STATES printf("%c: %s\n", character, __FUNCTION__);
#else
#define DEBUG_MSG_STATES
#endif

////////////////////////////////////////////////////////////////////////
// main ////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////

int main(int argc, char *argv[]) {
   // Check arguments
   if (argc <= 1) {
      printf("Usage: %s [filename]\n", argv[0]);
      exit(1);
   }

   // Create new cfg object
   struct cfg_t cfg;

   // Try to open file
   cfg.fhandle = fopen(argv[argc-1],"r");
   if (cfg.fhandle == NULL) {
      // Error
      printf("ERROR: Could not open file %s\n",argv[argc-1]);
      exit(1);
   }

   // Give name
   cfg.name = argv[argc-1];

   // Initialize
   cfg.stat_no = 1;
   cfg.stat_no_allocd = 4;
   cfg.stat_len        = (size_t*) malloc(cfg.stat_no_allocd*sizeof(size_t));
   cfg.stat_len_allocd = (size_t*) malloc(cfg.stat_no_allocd*sizeof(size_t));
   cfg.stat_len[0] = 0;
   cfg.stat_len_allocd[0] = 16;
   cfg.stat = (char**) malloc(cfg.stat_no_allocd*sizeof(char*));
   cfg.stat[0] = (char*) malloc(cfg.stat_len_allocd[0]*sizeof(char));
   cfg.dyn_no = 0;
   cfg.dyn_no_allocd = 0;
   cfg.dyn_alt_no = NULL;
   cfg.dyn_alt_no_allocd = NULL;
   cfg.dyn_len = NULL;
   cfg.dyn_len_allocd = NULL;
   cfg.dyn = NULL;
   cfg.tag = NULL;
   cfg.tag_len = NULL;
   cfg.tag_len_allocd = NULL;

   // Begin parsing
   in_static(&cfg);

   // Close file
   fclose(cfg.fhandle);

   #ifdef DEBUG
   // Print
   int i, j;
   size_t c;
   for (i=0; i<cfg.stat_no; i++) {
      printf("stat[%d]: ", i);
      for (c=0; c<cfg.stat_len[i]; c++) {
         printf("%c",cfg.stat[i][c]);
      }
      printf("\n");
   }
   for (i=0; i<cfg.dyn_no; i++) {
      printf("tag[%d]: ", i);
      for (c=0; c<cfg.tag_len[i]; c++) {
         printf("%c", cfg.tag[i][c]);
      }
      printf("\n");
      for (j=0; j<cfg.dyn_alt_no[i]; j++) {
         printf("dyn[%d][%d]: ",i,j);
         for (c=0; c<cfg.dyn_len[i][j]; c++) {
            printf("%c",cfg.dyn[i][j][c]);
         }
         printf("\n");
      }
   }
   #endif

   // Write all different files
   write_dyn_files(&cfg);

   // Free allocated space
   cleanup(&cfg);

   return 0;
}

////////////////////////////////////////////////////////////////////////
// write_dyn_files /////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////

void write_dyn_files(struct cfg_t* cfg) {
   // DEBUG
   int i, j;
   int overtake;
   size_t c;
   FILE* outfile;
   char* outname;

   // Calculate reference array
   calculate_reference_array(cfg);

   // Check if the number of alternatives in dynamic parts with the same tag is consistent
   check_tag_alt_no(cfg);

   #ifdef DEBUG
   printf("reference = ");
   for (i = 0; i < cfg->dyn_no; i++) printf("%d ",cfg->reference[i]);
   printf("\n");
   printf("reference_normalized = ");
   for (i = 0; i < cfg->dyn_no; i++) printf("%d ",cfg->reference_normalized[i]);
   printf("\n");
   printf("different_tag_no = %d\n",cfg->different_tag_no);
   printf("dyn_no = %d\n",cfg->dyn_no);
   #endif

   // Find out how many possibilities for different alternatives there are
   int possibilities = 1;
   for (i=0; i<cfg->dyn_no; i++) {
      // Only count dynamic parts with the same tag once
      if (cfg->reference[i] == i) {
         possibilities = possibilities * cfg->dyn_alt_no[i];
      }
   }

   // Create array that holds the index of the alternatives
   int* indices = malloc(cfg->different_tag_no*sizeof(int));
   for (i=0; i<cfg->different_tag_no; i++) indices[i] = 0;

   // Check for too many possibilities >= 10^7
   if (possibilities > 999999) {
      printf("Error: Possibilities exceed 10^7-1 = 9 999 999, aborting");
      exit(1);
   }

   // Loop over all possible files
   for (i=0; i<possibilities; i++) {

      #ifdef DEBUG
      // Print the indices
      printf("indices:");
      for (j=0; j<cfg->different_tag_no; j++) printf("  %d", indices[j]);
      printf("\n");
      #endif

      // Open output file
                       //                  "."  #  "." exp  \0
      outname = malloc((strlen(cfg->name) + 1 + 7 + 1 + 3 + 1)*sizeof(char));
      sprintf(outname, "%s.%07d.exp", cfg->name, i);
      outfile = fopen(outname,"w+");
      free(outname);

      // Print the files
      for (j=0; j<cfg->stat_no; j++) {
         // Print static part
         for (c=0; c<cfg->stat_len[j]; c++)
            fprintf(outfile, "%c", cfg->stat[j][c]);

         // Print dynamic part (if present)
         if (j < cfg->dyn_no) {
            for (c=0; c<cfg->dyn_len[j][indices[cfg->reference_normalized[j]]]; c++) {
               putc(cfg->dyn[j][indices[cfg->reference_normalized[j]]][c], outfile);
            }
         }
      }

      // Close file
      fclose(outfile);

      // Increase the index
      overtake = 1;
      for (j=0; j<cfg->dyn_no; j++) {
         // Only increase for new tags
         if (cfg->reference[j] == j) {
            if (overtake == 1) {
               indices[cfg->reference_normalized[j]] ++;
               // If the index is too big, set it back to zero and
               // increase the next index
               if (indices[cfg->reference_normalized[j]] >= cfg->dyn_alt_no[j]) {
                  indices[cfg->reference_normalized[j]] = 0;
                  overtake = 1;
               } else {
                  overtake = 0;
               }
            }
         }
      }
   }

   // Free indices
   free(indices);

}

////////////////////////////////////////////////////////////////////////
// cleanup /////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////

void cleanup(struct cfg_t* cfg) {
   // Free all allocated space
   int i, j;

   // Static text variables
   for (i=0; i<cfg->stat_no; i++) {
      free(cfg->stat[i]);
   }
   free(cfg->stat);
   free(cfg->stat_len);
   free(cfg->stat_len_allocd);

   // Dynamic text variables
   for (i=0; i<cfg->dyn_no; i++) {
      for (j=0; j<cfg->dyn_alt_no[i]; j++) {
         free(cfg->dyn[i][j]);
      }
      free(cfg->dyn[i]);
      free(cfg->dyn_len[i]);
      free(cfg->dyn_len_allocd[i]);
   }
   free(cfg->dyn);
   free(cfg->dyn_len);
   free(cfg->dyn_len_allocd);
   free(cfg->dyn_alt_no);
   free(cfg->dyn_alt_no_allocd);

   // Tag variables
   for (i=0; i<cfg->dyn_no; i++) {
      free(cfg->tag[i]);
   }
   free(cfg->tag);
   free(cfg->tag_len);
   free(cfg->tag_len_allocd);

   // Variables for reference
   free(cfg->reference);
   free(cfg->reference_normalized);

}

////////////////////////////////////////////////////////////////////////
// Helper functions ////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////

void check_dyn_len(struct cfg_t* cfg) {
   // Check length
   if (    cfg->dyn_len[cfg->dyn_no-1][cfg->dyn_alt_no[cfg->dyn_no-1]-1]
        >= cfg->dyn_len_allocd[cfg->dyn_no-1][cfg->dyn_alt_no[cfg->dyn_no-1]-1]) {
      // Reallocate if necessary
      if (cfg->dyn_len_allocd[cfg->dyn_no-1][cfg->dyn_alt_no[cfg->dyn_no-1]-1] == 0) {
         cfg->dyn_len_allocd[cfg->dyn_no-1][cfg->dyn_alt_no[cfg->dyn_no-1]-1] = 1;
      } else {
         cfg->dyn_len_allocd[cfg->dyn_no-1][cfg->dyn_alt_no[cfg->dyn_no-1]-1] *= 2;
      }
      cfg->dyn[cfg->dyn_no-1][cfg->dyn_alt_no[cfg->dyn_no-1]-1] =
         (char*) realloc(cfg->dyn[cfg->dyn_no-1][cfg->dyn_alt_no[cfg->dyn_no-1]-1],
                         cfg->dyn_len_allocd[cfg->dyn_no-1][cfg->dyn_alt_no[cfg->dyn_no-1]-1]*sizeof(char));
   }

   return;
}

void check_dyn_alt_no(struct cfg_t* cfg) {
   // Check length
   if (    cfg->dyn_alt_no[cfg->dyn_no-1]
        >= cfg->dyn_alt_no_allocd[cfg->dyn_no-1]) {
      // Reallocate if necessary
      if (cfg->dyn_alt_no_allocd[cfg->dyn_no-1] == 0) {
         cfg->dyn_alt_no_allocd[cfg->dyn_no-1] = 1;
      } else {
         cfg->dyn_alt_no_allocd[cfg->dyn_no-1] *= 2;
      }
      cfg->dyn[cfg->dyn_no-1] =
         (char**) realloc(cfg->dyn[cfg->dyn_no-1],
                          cfg->dyn_alt_no_allocd[cfg->dyn_no-1]*sizeof(char*));
      cfg->dyn_len[cfg->dyn_no-1] =
         (size_t*) realloc(cfg->dyn_len[cfg->dyn_no-1],
                           cfg->dyn_alt_no_allocd[cfg->dyn_no-1]*sizeof(size_t));
      cfg->dyn_len_allocd[cfg->dyn_no-1] =
         (size_t*) realloc(cfg->dyn_len_allocd[cfg->dyn_no-1],
                           cfg->dyn_alt_no_allocd[cfg->dyn_no-1]*sizeof(size_t));
   }

   return;
}

void check_dyn_no(struct cfg_t* cfg) {
   // Check length
   if (    cfg->dyn_no
        >= cfg->dyn_no_allocd) {
      // Reallocate if necessary
      if (cfg->dyn_no_allocd == 0) {
         cfg->dyn_no_allocd = 1;
      } else {
         cfg->dyn_no_allocd *= 2;
      }
      cfg->dyn_alt_no =
         (int*) realloc(cfg->dyn_alt_no,
                        cfg->dyn_no_allocd*sizeof(int));
      cfg->dyn_alt_no_allocd =
         (int*) realloc(cfg->dyn_alt_no_allocd,
                        cfg->dyn_no_allocd*sizeof(int));
      cfg->dyn_len =
         (size_t**) realloc(cfg->dyn_len,
                            cfg->dyn_no_allocd*sizeof(size_t*));
      cfg->dyn_len_allocd =
         (size_t**) realloc(cfg->dyn_len_allocd,
                            cfg->dyn_no_allocd*sizeof(size_t*));
      cfg->dyn =
         (char***) realloc(cfg->dyn,
                           cfg->dyn_no_allocd*sizeof(char**));
      cfg->tag =
         (char**) realloc(cfg->tag,
                          cfg->dyn_no_allocd*sizeof(char*));
      cfg->tag_len =
         (size_t*) realloc(cfg->tag_len,
                           cfg->dyn_no_allocd*sizeof(size_t));
      cfg->tag_len_allocd =
         (size_t*) realloc(cfg->tag_len_allocd,
                           cfg->dyn_no_allocd*sizeof(size_t));
   }

   return;
}

void alloc_next_dyn(struct cfg_t* cfg) {
          cfg->dyn_len[cfg->dyn_no-1][cfg->dyn_alt_no[cfg->dyn_no-1]-1] = 0;
   cfg->dyn_len_allocd[cfg->dyn_no-1][cfg->dyn_alt_no[cfg->dyn_no-1]-1] = 8;
   cfg->dyn[cfg->dyn_no-1][cfg->dyn_alt_no[cfg->dyn_no-1]-1] =
      (char*) malloc(cfg->dyn_len_allocd[cfg->dyn_no-1][cfg->dyn_alt_no[cfg->dyn_no-1]-1]*sizeof(char));
   return;
}

void alloc_next_dyn_alt(struct cfg_t* cfg) {
   // Dynamic text
   cfg->dyn_alt_no[cfg->dyn_no-1] = 0;
   cfg->dyn_alt_no_allocd[cfg->dyn_no-1] = 4;
   cfg->dyn[cfg->dyn_no-1] =
      (char**) malloc(cfg->dyn_alt_no_allocd[cfg->dyn_no-1]*sizeof(char*));
   cfg->dyn_len[cfg->dyn_no-1] =
      (size_t*) malloc(cfg->dyn_alt_no_allocd[cfg->dyn_no-1]*sizeof(size_t));
   cfg->dyn_len_allocd[cfg->dyn_no-1] =
      (size_t*) malloc(cfg->dyn_alt_no_allocd[cfg->dyn_no-1]*sizeof(size_t));

   // Tags
   cfg->tag_len[cfg->dyn_no-1] = 0;
   cfg->tag_len_allocd[cfg->dyn_no-1] = 8;
   cfg->tag[cfg->dyn_no-1] =
      (char*) malloc(cfg->tag_len_allocd[cfg->dyn_no-1]*sizeof(char));

   return;
}

void check_tag_len(struct cfg_t* cfg) {
   // Check length
   if (    cfg->tag_len[cfg->dyn_no-1]
        >= cfg->tag_len_allocd[cfg->dyn_no-1]) {
      // Reallocate if necessary
      cfg->tag_len_allocd[cfg->dyn_no-1] *= 2;
      cfg->tag[cfg->dyn_no-1] =
         (char*) realloc(cfg->tag[cfg->dyn_no-1],
                         cfg->tag_len_allocd[cfg->dyn_no-1]*sizeof(char));
   }

   return;
}

void check_stat_len(struct cfg_t* cfg) {
   // Check length
   if (     cfg->stat_len[cfg->stat_no-1]
         >= cfg->stat_len_allocd[cfg->stat_no-1]) {
      // Reallocate if necessary
      cfg->stat_len_allocd[cfg->stat_no-1] *= 2;
      cfg->stat[cfg->stat_no-1] =
         (char*) realloc(cfg->stat[cfg->stat_no-1],
                         cfg->stat_len_allocd[cfg->stat_no-1]*sizeof(char));
   }

   return;
}

void check_stat_no(struct cfg_t* cfg) {
   // Check number
   if ( cfg->stat_no >= cfg->stat_no_allocd) {
      // Reallocate if necessary
      cfg->stat_no_allocd *= 2;
      cfg->stat = (char**) realloc(cfg->stat, cfg->stat_no_allocd*sizeof(char*));
      cfg->stat_len = (size_t*) realloc(cfg->stat_len, cfg->stat_no_allocd*sizeof(size_t));
      cfg->stat_len_allocd = (size_t*) realloc(cfg->stat_len_allocd, cfg->stat_no_allocd*sizeof(size_t));
   }

   return;
}

void alloc_next_stat(struct cfg_t* cfg) {
   cfg->stat_len[cfg->stat_no-1] = 0;
   cfg->stat_len_allocd[cfg->stat_no-1] = 16;
   cfg->stat[cfg->stat_no-1] = (char*) malloc(cfg->stat_len_allocd[cfg->stat_no-1]*sizeof(char));
   return;
}

int read_character(struct cfg_t* cfg, char* character) {
   return fscanf(cfg->fhandle, "%c", character);
}

void read_character_err(struct cfg_t* cfg, char* character, char* err_msg) {
   if (read_character(cfg, character) < 1) {
      printf("ERROR: %s\n", err_msg);
      exit(1);
   }

   return;
}

#define read_character_suc(CFG, CHARACTER) \
if (read_character(CFG, CHARACTER) < 1) {  \
   return;                                  \
}

void write_static(struct cfg_t* cfg, char character) {
   // Check length of buffer and reallocate if it's too small
   check_stat_len(cfg);
   // Save character
   cfg->stat[cfg->stat_no-1][cfg->stat_len[cfg->stat_no-1]] = character;
   cfg->stat_len[cfg->stat_no-1]++;

   return;
}

void write_dynamic(struct cfg_t* cfg, char character) {
   // Check length of buffer and reallocate if it's too small
   check_dyn_len(cfg);
   // Save character
   cfg->dyn_len[cfg->dyn_no-1][cfg->dyn_alt_no[cfg->dyn_no-1]-1]++;
   cfg->dyn[cfg->dyn_no-1][cfg->dyn_alt_no[cfg->dyn_no-1]-1][cfg->dyn_len[cfg->dyn_no-1][cfg->dyn_alt_no[cfg->dyn_no-1]-1]-1]
      = character;

   return;
}

void write_tag(struct cfg_t* cfg, char character) {
   // Check length of buffer and reallocate if it's too small
   check_tag_len(cfg);
   // Save character
   cfg->tag[cfg->dyn_no-1][cfg->tag_len[cfg->dyn_no-1]] = character;
   cfg->tag_len[cfg->dyn_no-1]++;

   return;
}

void next_static_section(struct cfg_t* cfg) {
   // Check lenght of buffer and reallocate if it's too small
   check_stat_no(cfg);
   cfg->stat_no ++;
   // Allocate buffer
   alloc_next_stat(cfg);
}

void next_dynamic_alternative(struct cfg_t* cfg) {
   // Check lenght of buffer and reallocate if it's too small
   check_dyn_alt_no(cfg);
   cfg->dyn_alt_no[cfg->dyn_no-1]++;
   // Allocate buffer
   alloc_next_dyn(cfg);

   return;
}

void next_dynamic_section(struct cfg_t* cfg) {
   // Check lenght of buffer and reallocate if it's too small
   check_dyn_no(cfg);
   cfg->dyn_no++;
   // Allocate buffer
   alloc_next_dyn_alt(cfg);
   // Check lenght of buffer and reallocate if it's too small
   check_dyn_alt_no(cfg);
   cfg->dyn_alt_no[cfg->dyn_no-1]++;
   // Allocate buffer
   alloc_next_dyn(cfg);

   return;
}

int compare_tags(struct cfg_t* cfg, int index1, int index2) {
   // Check if one tag is empty (empty tags are not equal). If so return return 0 (false)
   if (   cfg->tag_len[index1] == 0
       || cfg->tag_len[index2] == 0) {
      return 0;
   }
   // Check if lengths are different. If so return 0 (false)
   if (cfg->tag_len[index1] != cfg->tag_len[index2]) {
      return 0;
   } else {
      size_t c;
      // Check all characters, if one is different, return 0 (false)
      for (c=0; c<cfg->tag_len[index1]; c++) {
         if (cfg->tag[index1][c] != cfg->tag[index2][c]) {
            return 0;
         }
      }
      // All characters are the same, so the tags are the same, return 1 (true)
      return 1;
   }
}

void calculate_reference_array(struct cfg_t* cfg) {
   // Allocate space for reference array
   cfg->reference = (int*) malloc(cfg->dyn_no*sizeof(int));
   cfg->reference_normalized = (int*) malloc(cfg->dyn_no*sizeof(int));

   // Count numbers that don't appear anymore
   int removed_no = 0;

   // Loop over all tags decreasing
   int i, j, remove_i;
   for (i = cfg->dyn_no-1; i > 0; i--) {
      // Default:
      cfg->reference[i] = i;
      cfg->reference_normalized[i] = i;

      // We need to know if an index was overwritten
      remove_i = 0;
      // Loop over all smaller indices, thus finding the foremost equal tag
      for (j = i-1; j >= 0; j--) {
         if (compare_tags(cfg, i, j)) {
            // Tags with indices i and j are the same and j<=i
            cfg->reference[i] = j;
            cfg->reference_normalized[i] = j;
            remove_i = 1;
         }
      }

      // If index was overwritten,
      if (remove_i) {
         // Increase number of removed indices
         removed_no++;


         // Decrease all following indices greater than i for the normalized reference
         for (j = i+1; j < cfg->dyn_no; j++) {
            if (cfg->reference_normalized[j]>i) cfg->reference_normalized[j]--;
         }

      }
   }
   // The first index has to stay
   cfg->reference[0] = 0;
   cfg->reference_normalized[0] = 0;


   // Calculate number of different tags
   cfg->different_tag_no = cfg->dyn_no - removed_no;

   return;
}

void check_tag_alt_no(struct cfg_t* cfg) {
   // Loop over all dynamic parts
   int i;
   for (i = 0; i < cfg->dyn_no; i++) {
      // Check if there is an alternative with the same tag
      if (cfg->reference[i] != i) {
         // Check if the number of alternatives is the same
         if (cfg->dyn_alt_no[i] != cfg->dyn_alt_no[cfg->reference[i]]) {
            // If not, error
            printf("Error: Number of alternatives in the dynamic parts with a common tag is not consistent.\n");
            exit(1);
         }
      }
   }

   return;
}

////////////////////////////////////////////////////////////////////////
// Begin of the state functions ////////////////////////////////////////
////////////////////////////////////////////////////////////////////////

void in_static(struct cfg_t* cfg) {
   // Read one character from the file (success if EOF)
   char character;
   read_character_suc(cfg, &character);

   DEBUG_MSG_STATES

   // Test the character
   if (character == '[') {
      have_opening_bracket(cfg);
   }
   else {
      // Character is static
      write_static(cfg, character);

      in_static(cfg);
   }

   return;
}

void have_opening_bracket(struct cfg_t* cfg) {
   // Read one character from the file (success if EOF)
   char character;
   read_character_suc(cfg, &character)

   DEBUG_MSG_STATES

   // Test the character
   if (character == '[') {
      // We have found a new dynamic part
      next_dynamic_section(cfg);

      first_in_dynamic(cfg);
   }
   else {
      // Characters '[' and character are static
      write_static(cfg, '[');
      write_static(cfg, character);

      in_static(cfg);
   }

   return;
}

void first_in_dynamic(struct cfg_t* cfg) {
   // Read one character from the file (Error if EOF)
   char character;
   read_character_err(cfg, &character, "Unmatched [[");

   DEBUG_MSG_STATES

   // Test the character
   if (character == '(') {
      have_opening_paren(cfg);
   } else {
      // Character is dynamic
      write_dynamic(cfg, character);

      in_dynamic(cfg);
   }

   return;
}

void have_opening_paren(struct cfg_t* cfg) {
   // Read one character from the file (Error if EOF)
   char character;
   read_character_err(cfg, &character, "Unmatched [[");

   DEBUG_MSG_STATES

   // Test the character
   if (character == '(') {
      // TODO

      in_tag(cfg);
   } else {
      // Character '(' was dynamic
      write_dynamic(cfg, '(');
      // Character character is dynamic
      write_dynamic(cfg, character);

      in_dynamic(cfg);
   }

   return;
}

void in_tag(struct cfg_t* cfg) {
   // Read one character from the file (Error if EOF)
   char character;
   read_character_err(cfg, &character, "Unmatched (( and therefore unmatched [[, as well");

   DEBUG_MSG_STATES

   // Test the character
   if (character == ')') {
      have_closing_paren(cfg);
   } else {
      // Character is part of the tag
      write_tag(cfg, character);

      in_tag(cfg);
   }

   return;
}

void have_closing_paren(struct cfg_t* cfg) {
   // Read one character from the file (Error if EOF)
   char character;
   read_character_err(cfg, &character, "Unmatched (( and therefore umatched [[, as well");

   DEBUG_MSG_STATES

   // Test the character
   if (character == ')') {
      in_dynamic(cfg);
   } else {
      // Character ')' was part of the tag
      write_tag(cfg, ')');
      // Character character is part of the tag
      write_tag(cfg, character);

      in_tag(cfg);
   }

   return;
}

void in_dynamic(struct cfg_t* cfg) {
   // Read one character from the file (Error if EOF)
   char character;
   read_character_err(cfg, &character, "Unmatched [[");

   DEBUG_MSG_STATES

   // Test the character
   if (character == '|') {
      have_delimiter(cfg);
   }
   else if (character == ']') {
     have_closing_bracket(cfg);
   }
   else {
      // Character is dynamic
      write_dynamic(cfg, character);

      in_dynamic(cfg);
   }

   return;
}


void have_delimiter(struct cfg_t* cfg) {
   // Read one character from the file (Error if EOF)
   char character;
   read_character_err(cfg, &character, "Unmatched [[");

   DEBUG_MSG_STATES

   // Test the character
   if (character == '|') {
      // Next alternative
      next_dynamic_alternative(cfg);

      in_dynamic(cfg);
   }
   else {
      // Character '|' was normal text
      write_dynamic(cfg, '|');
      // Character character is dynamic
      write_dynamic(cfg, character);

      in_dynamic(cfg);
   }

   return;
}

void have_closing_bracket(struct cfg_t* cfg) {
   // Read one character from the file (Error if EOF)
   char character;
   read_character_err(cfg, &character, "Unmatched ]]");

   DEBUG_MSG_STATES

   // Test the character
   if (character == ']') {
      // New static section
      next_static_section(cfg);

      in_static(cfg);
   }
   else {
      // Character ']' was normal text
      write_dynamic(cfg, ']');
      // Character character is dynamic
      write_dynamic(cfg, character);

      in_dynamic(cfg);
   }

   return;
}
