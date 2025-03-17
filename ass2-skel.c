/* Program to perform 1D cellular automaton (CA) computations and to use 1D CA
   to solve the density classification problem.

  Skeleton program written by Artem Polyvyanyy, http://polyvyanyy.com/,
  September 2024, with the intention that it be modified by students
  to add functionality, as required by the assignment specification.
  All included code is (c) Copyright University of Melbourne, 2024.

  Authorship Declaration:

  (1) I certify that except for the code provided in the initial skeleton file,
  the program contained in this submission is completely my own individual
  work, except where explicitly noted by further comments that provide details
  otherwise. I understand that work that has been developed by another student,
  or by me in collaboration with other students, or by non-students as a result
  of request, solicitation, or payment, may not be submitted for assessment in
  this subject. I understand that submitting for assessment work developed by
  or in collaboration with other students or non-students constitutes Academic
  Misconduct, and may be penalized by mark deductions, or by other penalties
  determined via the University of Melbourne Academic Honesty Policy, as
  described at https://academicintegrity.unimelb.edu.au.

  (2) I also certify that I have not provided a copy of this work in either
  softcopy or hardcopy or any other form to any other student, and nor will I
  do so until after the marks are released. I understand that providing my work
  to other students, regardless of my intention or any undertakings made to me
  by that other student, is also Academic Misconduct.

  (3) I further understand that providing a copy of the assignment specification
  to any form of code authoring or assignment tutoring service, or drawing the
  attention of others to such services and code that may have been made
  available via such a service, may be regarded as Student General Misconduct
  (interfering with the teaching activities of the University and/or inciting
  others to commit Academic Misconduct). I understand that an allegation of
  Student General Misconduct may arise regardless of whether or not I personally
  make use of such solutions or sought benefit from such actions.

  Signed by: Eric Zhang
  Dated:     10/10/24
*/

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>

/* #DEFINE'S -----------------------------------------------------------------*/
#define SDELIM "==STAGE %d============================\n" // stage delimiter
#define MDELIM "-------------------------------------\n"  // delimiter of -'s
#define THEEND "==THE END============================\n"  // end message

#define CRTRNC '\r' // carriage return character
#define NBRHDS 8    // number of possible neighborhoods
#define PATTERN 3

#define STAGE_0 0
#define STAGE_1 1
#define STAGE_2 2
#define INITIAL 0

#define ON '1'
#define OFF '0'
#define INT_CONVERSION '0'

#define FIRST_BIT 1
#define SECOND_BIT 2
#define THIRD_BIT 4

#define RULE_184 ((unsigned char)184)
#define RULE_232 ((unsigned char)232)

/* TYPE DEFINITIONS ----------------------------------------------------------*/
typedef char cells_t;                 // base type to store states of cells
typedef struct state state_t;         // a cellular automaton state
typedef unsigned char rule_t[NBRHDS]; // an elementary CA update rule function

typedef struct
{
    int on_count;
    int off_count;
    int cell_no;
    int start_time;
} run_data_t;

struct state
{                  // a state in a CA is defined by
    cells_t *clls; // ... an array of cells and
    state_t *next; // ... a link to the next state
};

typedef struct
{                  // a run of a CA consists of
    state_t *init; // ... the initial state and
    state_t *curr; // ... the current state,
} run_t;           // implemented as a linked list of states

typedef struct
{                       // an elementary CA is defined by
    unsigned char code; // ... a code of the update rule,
    unsigned int size;  // ... a number of cells,
    unsigned int time;  // ... the current time step,
    rule_t rule;        // ... an update rule function, and
    run_t *run;         // ... a run of state steps
} CA_t;

/* USEFUL FUNCTIONS ----------------------------------------------------------*/
int mygetchar(void); // getchar() that skips
                     //    carriage returns
CA_t *read_initial_configuration(void);
CA_t *initialise_ca(void);
void read_initial_cell_state(CA_t *ca);
void convert_code_to_rule_t(CA_t *ca, unsigned char code);
unsigned char *decimal_to_binary(int decimal, int num_bits);
void print_stage_zero(CA_t *ca);
void print_neighbourhoods(void);
void print_cell_state(cells_t *clls, int size, int state_no);
int get_rule_index(int i, CA_t *ca);
cells_t *evolve_cell_step(CA_t *ca);
run_t *insert_next_evolution_step(run_t *run, cells_t *evolved_cells);
void execute_automata(CA_t *ca, int time_steps);
run_data_t *initialise_run_data(void);
state_t *find_start_time(CA_t *ca, int start_time);
void count_on_off(run_data_t *run_data, state_t *run_ptr);
void print_run_data(run_data_t *run_data, CA_t *ca);
void print_stage_one(run_data_t **run_data, CA_t *ca);
CA_t *create_next_automata(CA_t *curr_ca, unsigned char code);
void on_off_at_time_step(CA_t *ca);
void print_stage_two(CA_t **ca184, CA_t **ca232, CA_t *ca,
                     run_data_t **stage_two_run);
void free_run(run_t *run);
void free_memory(CA_t **ca, CA_t **ca184, CA_t **ca232,
                 run_data_t **stage_one_run, run_data_t **stage_two_run);

/* WHERE IT ALL HAPPENS ------------------------------------------------------*/
int main(int argc, char *argv[])
{

    CA_t *ca = read_initial_configuration();
    print_stage_zero(ca);

    run_data_t *stage_one_run = NULL;
    print_stage_one(&stage_one_run, ca);

    CA_t *ca184 = NULL;
    CA_t *ca232 = NULL;
    run_data_t *stage_two_run = NULL;
    print_stage_two(&ca184, &ca232, ca, &stage_two_run);
    free_memory(&ca, &ca184, &ca232, &stage_one_run, &stage_two_run);

    return EXIT_SUCCESS; // algorithms are fun!!!
}

/* USEFUL FUNCTIONS ----------------------------------------------------------*/

// An improved version of getchar(); skips carriage return characters.
// NB: Adapted version of the mygetchar() function by Alistair Moffat
int mygetchar()
{

    int c;
    while ((c = getchar()) == CRTRNC)
        ; // skip carriage return characters
    return c;
}

/* THE END -------------------------------------------------------------------*/

// Reads from input and configures components of initial cellular automata
CA_t *read_initial_configuration(void)
{

    CA_t *ca;
    ca = initialise_ca();

    // Read size and code from input
    scanf("%u\n%hhu\n", &ca->size, &ca->code);
    // Set initial time to 0
    ca->time = 0;

    // Allocate memory for cell array
    ca->run->init->clls = (cells_t *)malloc((ca->size) * sizeof(cells_t));
    assert(ca->run->init->clls != NULL);

    read_initial_cell_state(ca);

    convert_code_to_rule_t(ca, ca->code);

    return ca;
}

/*-------------------------------------------------------------------*/

// Allocates memory for cellular automata structure
CA_t *initialise_ca(void)
{
    // Malloc ca
    CA_t *ca = (CA_t *)malloc(sizeof(CA_t));
    assert(ca != NULL);

    // Malloc run
    ca->run = (run_t *)malloc(sizeof(run_t));
    assert(ca->run != NULL);

    // Malloc state
    ca->run->init = (state_t *)malloc(sizeof(state_t));
    assert(ca->run->init != NULL);

    // link current and init states
    ca->run->curr = ca->run->init;

    return ca;
}

/*-------------------------------------------------------------------*/

// Reads initial cell state from input and stores in cells component of ca
void read_initial_cell_state(CA_t *ca)
{

    for (int i = 0; i < (ca->size); i++)
    {
        char ch;
        ch = mygetchar();

        if (ch == '*')
        {
            // define magic number
            ca->run->init->clls[i] = ON;
        }
        else if (ch == '.')
        {
            // define magic number
            ca->run->init->clls[i] = OFF;
        }
    }
}

/*-------------------------------------------------------------------*/

// Converts code to binary string representation and stores in rule
void convert_code_to_rule_t(CA_t *ca, unsigned char code)
{

    unsigned char *binary = decimal_to_binary(code, NBRHDS);

    for (int i = 0; i < NBRHDS; i++)
    {
        ca->rule[NBRHDS - 1 - i] = binary[i];
    }
    free(binary);
}

/*-------------------------------------------------------------------*/

// Converts integer decimal into a binary string representation
unsigned char *decimal_to_binary(int decimal, int num_bits)
{

    // Allocates memory array for binary string representation
    unsigned char *binary = (unsigned char *) \
                            malloc(num_bits * sizeof(unsigned char *));

    assert(binary != NULL);

    // Converts decimal to binary and populates binary array
    for (int i = num_bits - 1; i >= 0; i--)
    {
        if (decimal & (1 << i))
        {
            binary[num_bits - 1 - i] = ON;
        }
        else
        {
            binary[num_bits - 1 - i] = OFF;
        }
    }
    binary[num_bits] = '\0';

    return binary;
}

/*-------------------------------------------------------------------*/

// Print required output for stage zero
void print_stage_zero(CA_t *ca)
{

    printf(SDELIM, STAGE_0);

    printf("SIZE: %u\n", ca->size);
    printf("RULE: %hhu\n", ca->code);

    printf(MDELIM);
    print_neighbourhoods();

    // Prints rule for stage zero
    for (int j = 0; j < NBRHDS; j++)
    {
        printf("  %c ", ca->rule[j]);
    }
    printf("\n");

    printf(MDELIM);

    // Print initial cell state
    print_cell_state(ca->run->curr->clls, ca->size, ca->time);
}

/*-------------------------------------------------------------------*/

// Prints neighbourhood patterns for stage zero
void print_neighbourhoods()
{

    for (int i = 0; i < NBRHDS; i++)
    {
        unsigned char *binary = decimal_to_binary(i, PATTERN);
        printf(" %s", binary);
        free(binary);
    }
    printf("\n");
}

/*-------------------------------------------------------------------*/

// Prints time step and corresponding cell array
void print_cell_state(cells_t *clls, int size, int time_step)
{

    printf("%4d: ", time_step);

    for (int k = 0; k < size; k++)
    {
        if (clls[k] == ON)
        {
            printf("*");
        }
        else if (clls[k] == OFF)
        {
            printf(".");
        }
    }
    printf("\n");
}

/*-------------------------------------------------------------------*/
// Gets index for rule that matches with the neighbourhood pattern
int get_rule_index(int i, CA_t *ca)
{

    // Find left middle and right cells
    // and converts characeters to respective integers
    int left = ca->run->curr->clls[(i - 1 + ca->size) % ca->size] \
               - INT_CONVERSION;

    int middle = ca->run->curr->clls[i] - INT_CONVERSION;
    int right = ca->run->curr->clls[(i + 1) % ca->size] - INT_CONVERSION;

    // Converts binary pattern to decimal integer
    return (left * THIRD_BIT + middle * SECOND_BIT + right * FIRST_BIT);
}

/*-------------------------------------------------------------------*/
// Returns evolved cell array based on rule
cells_t *evolve_cell_step(CA_t *ca)
{

    // Allocate memory for evolved cell array
    cells_t *evolved_cells = (cells_t *)malloc((ca->size) * sizeof(cells_t *));
    assert(evolved_cells != NULL);

    for (int i = 0; i < (ca->size); i++)
    {
        evolved_cells[i] = ca->rule[get_rule_index(i, ca)];
    }

    return evolved_cells;
}

/*-------------------------------------------------------------------*/

// Inserts the next evolution step of cells at the end of run
// NB: Adapted version of the insert_at_foot() function by Alistair Moffat
run_t *
insert_next_evolution_step(run_t *run, cells_t *evolved_cells)
{
    // Malloc new space for next_step
    state_t *next_step = (state_t *)malloc(sizeof(*next_step));
    assert(run != NULL && next_step != NULL);

    // Insert data into insert_next_step
    next_step->clls = evolved_cells;
    next_step->next = NULL;

    // Point curr to next_step
    run->curr->next = next_step;
    run->curr = next_step;

    return run;
}

/*-------------------------------------------------------------------*/

// Executes automata and evolves for specified time steps
void execute_automata(CA_t *ca, int time_steps)
{

    // Print last cell array from previous stage
    print_cell_state(ca->run->curr->clls, ca->size, ca->time);

    // Prints cell array for each evolution step
    while (ca->time < time_steps)
    {
        ca->time += 1;
        insert_next_evolution_step(ca->run, evolve_cell_step(ca));
        print_cell_state(ca->run->curr->clls, ca->size, ca->time);
    }
}

/*-------------------------------------------------------------------*/

// Allocates memory for run data structure and reads from input to
// configure run start time and cell number in array
run_data_t *initialise_run_data(void)
{

    run_data_t *run_data = (run_data_t *)malloc(sizeof(run_data_t));
    run_data->on_count = 0;
    run_data->off_count = 0;

    scanf("%d,%d\n", &run_data->cell_no, &run_data->start_time);

    return run_data;
}

/*-------------------------------------------------------------------*/

// Starts from init of run and traverses list until start time in run is found
state_t *find_start_time(CA_t *ca, int start_time)
{

    // Initialise a pointer and counter to traverse through linked list
    state_t *run_ptr = ca->run->init;
    int stage_counter = 0;

    while (run_ptr != NULL && stage_counter < start_time)
    {
        run_ptr = run_ptr->next;
        stage_counter++;
    }
    return run_ptr;
}

/*-------------------------------------------------------------------*/

// Count number of ON and OFF states for run and stores in run data struct
void count_on_off(run_data_t *run_data, state_t *run_ptr)
{

    while (run_ptr != NULL)
    {
        if (run_ptr->clls[run_data->cell_no] == ON)
        {
            run_data->on_count += 1;
        }
        else if (run_ptr->clls[run_data->cell_no] == OFF)
        {
            run_data->off_count += 1;
        }
        run_ptr = run_ptr->next;
    }
}

/*-------------------------------------------------------------------*/

// Prints on and off count for cell number in array from a
// specified start time for a given run
void print_run_data(run_data_t *run_data, CA_t *ca)
{

    printf(MDELIM);
    count_on_off(run_data, find_start_time(ca, run_data->start_time));
    printf("#ON=%d #OFF=%d CELL#%d START@%d\n", run_data->on_count,
           run_data->off_count, run_data->cell_no, run_data->start_time);
}

/*-------------------------------------------------------------------*/

// Prints required output for stage one
void print_stage_one(run_data_t **run_data, CA_t *ca)
{

    printf(SDELIM, STAGE_1);

    int time_steps;
    scanf("%d\n", &time_steps);

    *run_data = initialise_run_data();
    execute_automata(ca, time_steps);
    print_run_data(*run_data, ca);
}

/*-------------------------------------------------------------------*/

// Allocates memory and configures a cellular automata structure
// based on the preceding automata. Links multiple automata with existing run
CA_t *create_next_automata(CA_t *curr_ca, unsigned char new_code)
{
    // Allocate memory for cellular automata structure
    CA_t *new_ca = (CA_t *)malloc(sizeof(CA_t));
    assert(new_ca != NULL);
    // malloc state

    // Initialise components
    new_ca->size = curr_ca->size;
    new_ca->code = new_code;
    new_ca->time = curr_ca->time;

    convert_code_to_rule_t(new_ca, new_code);

    // Reuse the existing run structure
    new_ca->run = curr_ca->run;

    return new_ca;
}

/*-------------------------------------------------------------------*/

// Count on and off cells for one cell array at a given time step
void on_off_at_time_step(CA_t *ca)
{

    cells_t *time_step_cells = find_start_time(ca, ca->time)->clls;
    int on_counter = 0;
    int off_counter = 0;

    for (int i = 0; i < ca->size; i++)
    {
        char cell = time_step_cells[i];
        if (cell == ON)
        {
            on_counter += 1;
        }
        else if (cell == OFF)
        {
            off_counter += 1;
        }
    }

    print_cell_state(time_step_cells, ca->size, ca->time);

    // Prints output based on on and off counters
    if (on_counter > off_counter)
    {
        printf("AT T=%d: #ON/#CELLS > 1/2\n", ca->time);
    }
    else if (on_counter < off_counter)
    {
        printf("AT T=%d: #ON/#CELLS < 1/2\n", ca->time);
    }
    else
    {
        printf("AT T=%d: #ON/#CELLS = 1/2\n", ca->time);
    }
}

/*-------------------------------------------------------------------*/

// Print required output for stage two
void print_stage_two(CA_t **ca184, CA_t **ca232, CA_t *ca,
                     run_data_t **stage_two_run)
{

    // Calculates time steps for stage 2
    int n, m;
    n = (ca->size - 2) / 2;
    m = (ca->size - 1) / 2;

    printf(SDELIM, STAGE_2);

    printf("RULE: %d; STEPS: %d.\n", RULE_184, n);
    printf(MDELIM);

    // Execute cellular automata with rule 184
    *ca184 = create_next_automata(ca, RULE_184);
    execute_automata(*ca184, ca->time + n);

    printf(MDELIM);

    printf("RULE: %d; STEPS: %d.\n", RULE_232, m);
    printf(MDELIM);

    // Execute cellular automata with rule 232
    *ca232 = create_next_automata(*ca184, RULE_232);
    execute_automata(*ca232, (*ca184)->time + m);

    *stage_two_run = initialise_run_data();
    print_run_data(*stage_two_run, ca);
    printf(MDELIM);

    on_off_at_time_step(ca);

    printf(THEEND);
}

/*-------------------------------------------------------------------*/

// Free memory for run component in cellular automata structures
// NB: Adapted version of the free_list() function by Alistair Moffat
void free_run(run_t *run)
{

    assert(run != NULL);
    // Free states in run
    state_t *curr = run->init, *prev;
    while (curr)
    {
        prev = curr;
        curr = curr->next;
        free(prev->clls);
        free(prev);
    }
    // Free the run itself
    free(run);
}

/*-------------------------------------------------------------------*/

// Frees all memory structures allocated in program
void free_memory(CA_t **ca, CA_t **ca184, CA_t **ca232,
                 run_data_t **stage_one_run, run_data_t **stage_two_run)
{

    // Free run_data structs
    assert(*stage_one_run != NULL);
    free(*stage_one_run);

    assert(*stage_two_run != NULL);
    free(*stage_two_run);

    // Free the original cellular automata structure and its run
    assert(*ca != NULL);
    assert((*ca)->run != NULL);
    free_run((*ca)->run);
    free(*ca);

    // Free remaining cellular automata structures
    assert(*ca184 != NULL);
    free(*ca184);
    assert(*ca232 != NULL);
    free(*ca232);

    // Sets all deallocated pointers to null
    *ca = NULL;
    *ca184 = NULL;
    *ca232 = NULL;
    *stage_one_run = NULL;
    *stage_two_run = NULL;
}
