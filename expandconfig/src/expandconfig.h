#include <stdio.h>
#include <stdlib.h>
#include <string.h>

struct cfg_t {
   // File handle
   FILE* fhandle;
   // Name
   char* name;
   // Number of static parts
   int stat_no;
   int stat_no_allocd;
   // Lengths of static parts in the file
   size_t* stat_len;
   size_t* stat_len_allocd;
   // Static parts
   char**  stat;
   // Number of dynamic parts
   int dyn_no;
   int dyn_no_allocd;
   // Number of alternatives in dynamic part
   int* dyn_alt_no;
   int* dyn_alt_no_allocd;
   // Length of alternatives of dynamic parts
   size_t** dyn_len;
   size_t** dyn_len_allocd;
   // Dynamic parts
   char***  dyn;
   // Tags of dynamic parts
   char** tag;
   // Length of the tags
   size_t* tag_len;
   size_t* tag_len_allocd;
   // Reference array
   int* reference;
   int* reference_normalized;
   // Number of different tags
   int different_tag_no;
};

int main(int argc, char *argv[]);
void write_dyn_files(struct cfg_t* cfg);
void cleanup(struct cfg_t* cfg);
void check_dyn_len(struct cfg_t* cfg);
void check_dyn_alt_no(struct cfg_t* cfg);
void check_dyn_no(struct cfg_t* cfg);
void alloc_next_dyn(struct cfg_t* cfg);
void alloc_next_dyn_alt(struct cfg_t* cfg);
void check_tag_len(struct cfg_t* cfg);
void check_stat_len(struct cfg_t* cfg);
void check_stat_no(struct cfg_t* cfg);
void alloc_next_stat(struct cfg_t* cfg);
int read_character(struct cfg_t* cfg, char* character);
void read_character_err(struct cfg_t* cfg, char* character, char* err_msg);
void write_static(struct cfg_t* cfg, char character);
void write_dynamic(struct cfg_t* cfg, char character);
void next_static_section(struct cfg_t* cfg);
void next_dynamic_alternative(struct cfg_t* cfg);
void next_dynamic_section(struct cfg_t* cfg);
int compare_tags(struct cfg_t* cfg, int index1, int index2);
void calculate_reference_array(struct cfg_t* cfg);
void check_tag_alt_no(struct cfg_t* cfg);
void in_static(struct cfg_t* cfg);
void have_opening_bracket(struct cfg_t* cfg);
void first_in_dynamic(struct cfg_t* cfg);
void have_opening_paren(struct cfg_t* cfg);
void in_tag(struct cfg_t* cfg);
void have_closing_paren(struct cfg_t* cfg);
void in_dynamic(struct cfg_t* cfg);
void have_delimiter(struct cfg_t* cfg);
void have_closing_bracket(struct cfg_t* cfg);
