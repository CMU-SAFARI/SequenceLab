/*
 * Copyright (c) <2016 - 2020>, Bilkent University and ETH Zurich 
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without modification,
 * are permitted provided that the following conditions are met:
 *
 * Redistributions of source code must retain the above copyright notice, this list
 * of conditions and the following disclaimer.
 * - Redistributions in binary form must reproduce the above copyright notice, this
 *   list of conditions and the following disclaimer in the documentation and/or other
 *   materials provided with the distribution.
 * - Neither the names of the Bilkent University, ETH Zurich,
 *   nor the names of its contributors may be
 *   used to endorse or promote products derived from this software without specific
 *   prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 * A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR
 * CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 * EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
 * PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
 * PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
 * LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 * NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

  Authors: 
  Mohammed Alser
	  mohammedalser AT bilkent DOT edu DOT tr
  Date:
  December 3rd, 2016
*/

#define GRID_LINE_LENGTH 512
#define GRID_LINES 1024
#define FILTER_COUNT 13
// if a filter takes more than 24 hours in total, we deactivate it.
#define TIME_LIMIT_IN_SECONDS 60*60*24

#define GRID_CSV_SUFFIX "edlib_estimate,error_threshold,read_length"

#define ADJACENCY_INDEX 0
#define BASE_COUNTING_INDEX 1
#define SHOUJI_INDEX 2
#define QGRAM_INDEX 3
#define MAGNET_INDEX 4
#define HD_INDEX 5
#define SHD_INDEX 6
#define SNEAKYSNAKE_INDEX 7
#define GRIM_INDEX 8
#define KSW2_INDEX 9
#define EDLIB_INDEX 10
#define GRIM_ORIGINAL_INDEX 11
#define QGRAM_LONG_INDEX 12


/* Include Files */
#include "main.h"
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <getopt.h>
#include <sys/stat.h>

#include "src/edlib.h"
#include <stdlib.h>
#include "src/needleman_wunsch.h"
#include <immintrin.h>
#include "filters/SneakySnake/SneakySnake.h" //SneakySnake for CPU
#include "filters/base-counting/Base_Counting.h" //mrFAST base counting
#include "filters/qgram/qgram.h"
#include "filters/shouji/Shouji.h"
#include "filters/hamming-distance/HD.h"
#include "filters/grim/grim.h"
#include "filters/pigeonhole/pigeonhole.h"
#include "aligners/ksw2/ksw2.h"
//#include "filters/swift/swift.h"
//#include "filters/GASSST/withgap.h" //GASSST filters (Tiled NW + Base Counting)
//#include "parasail.h"
//#include "parasail/matrices/blosum62.h"

/* to compile type the following 
sudo ldconfig -v
gcc -mavx2 -g -O3 -Wall -o main *.c -lz -lm -L libs/string_buffer/ -L src/ -ledlib -lalign -lstrbuf -lpthread
./main 0 4 8 100 /home/alser-xilinx/Desktop/Filters_29_11_2016/HumanGenome_Mason_mrfast_E15_300bp.txt 30000000 3
OR: use the following to check the memory leaks
valgrind --leak-check=yes --show-leak-kinds=all ./main
*/

static int full_edlib = 0;

static struct option long_options[] = {
        {"debug-mode", required_argument, NULL, 'd'},
        {"shouji-grid-size", required_argument, NULL, 's'},
        {"kmer-length", required_argument, NULL, 'k'},
        {"maxlength", required_argument, NULL, 'm'},
        {"file", required_argument, NULL, 'f'},
        {"read-count", required_argument, NULL, 'c'},
        {"iteration-count", required_argument, NULL, 'i'},
        {"grid-file-prefix", required_argument, NULL, 'g'},
        {"histogram-file-prefix", required_argument, NULL, 'h'},
        {"runtime-filename", required_argument, NULL, 'r'},
        {"max-edit-distance", required_argument, NULL, 'u'},
        {"min-edit-distance", required_argument, NULL, 'l'},
        {"activate-individually", required_argument, NULL, 'a'},
        {"full-edlib", no_argument, &full_edlib, 1},
        {"help", no_argument, NULL, 256},
        {0, 0, 0, 0}
};

char * linear_strcat(char * dest, char * src) {
    while (*dest) dest++;
    while ((*dest++ = *src++));
    return --dest;
}

void align(const char * tseq, const char * qseq, int sc_mch, int sc_mis, int gapo, int gape) {
    int i, a = sc_mch, b = sc_mis < 0? sc_mis : -sc_mis; // a>0 and b<0
    int8_t mat[25] = { a,b,b,b,0, b,a,b,b,0, b,b,a,b,0, b,b,b,a,0, 0,0,0,0,0 };
    int tl = strlen(tseq), ql = strlen(qseq);
    uint8_t *ts, *qs, c[256];
    ksw_extz_t ez;

    memset(&ez, 0, sizeof(ksw_extz_t));
    memset(c, 4, 256);
    c['A'] = c['a'] = 0; c['C'] = c['c'] = 1;
    c['G'] = c['g'] = 2; c['T'] = c['t'] = 3; // build the encoding table
    ts = (uint8_t*)malloc(tl);
    qs = (uint8_t*)malloc(ql);
    for (i = 0; i < tl; ++i) ts[i] = c[(uint8_t)tseq[i]]; // encode to 0/1/2/3
    for (i = 0; i < ql; ++i) qs[i] = c[(uint8_t)qseq[i]];
    ksw_extz(0, ql, qs, tl, ts, 5, mat, gapo, gape, -1, -1, 0, &ez);
//    for (i = 0; i < ez.n_cigar; ++i) // print CIGAR
//        printf("%d%c", ez.cigar[i]>>4, "MID"[ez.cigar[i]&0xf]);
//    putchar('\n');
    free(ez.cigar); free(ts); free(qs);
}

void csv_header(char * header, char ** filter_names, int use_trailing_comma) {
    int offset = 0;
    for (int i  = 0 ; i < FILTER_COUNT-1; i++) {
        sprintf(&header[offset], "%s,", filter_names[i]);
        offset += strlen(filter_names[i])+1; //+1 to account for comma
    }
    if (use_trailing_comma) {
        sprintf(&header[offset], "%s,", filter_names[FILTER_COUNT-1]);
    } else {
        sprintf(&header[offset], "%s", filter_names[FILTER_COUNT-1]);
    }
}

void print_help(int debug_mode, int shouji_grid_size, int adjacency_kmer_length, int read_length, char * tss_filename, int read_count, char * grid_filename, char * histogram_filename, char * runtime_filename, int max_edit_distance, int min_edit_distance) {
    printf("Usage: ./main [options]\n");
    printf("Required:\n");
    printf("\t-f, --file <filename>                     path to the .tss file to be analyzed\n");
    printf("Options:\n");
    printf("\t-d, --debug-mode <0|1>                    enable debug logs for each filter [%d]\n", debug_mode);
    printf("\t-s, --shouji-grid-size <int>              configure the grid size for the Shouji algorithm [%d]\n", shouji_grid_size);
    printf("\t-k, --adjacency-kmer-length <int>         k-mer length for the Adjacency filter [%d]\n", adjacency_kmer_length);
    printf("\t-m, --maxlength <int>                     truncate the reads to the given number of bases [%d]\n", read_length);
    printf("\t-c, --read-count <int>                    evaluate filters using only the given number of sequence pairs from the input file [%d]\n", read_count);
    printf("\t-g, --grid-file-prefix <filename>         use the given prefix for the grid files [\"%s\"]\n", grid_filename);
    printf("\t-h, --histogram-file-prefix <filename>    use the given prefix for the histogram file [\"%s\"]\n", histogram_filename);
    printf("\t-r, --runtime-filename <filename>         use the given prefix for the runtime file [\"%s\"]\n", runtime_filename);
    printf("\t-u, --max-edit-distance <int>             sweep the edit distance threshold up to the given maximum (in percent) [%d]\n", max_edit_distance);
    printf("\t-l, --min-edit-distance <int>             sweep the edit distance threshold from the given maximum (in percent) [%d]\n", min_edit_distance);
    printf("\t-a, --activate-individually <binary mask> evaluate only filters enabled in the given binary mask [1111111111111]\n");
    printf("\t                                          The mask applies from left to right to each filter in the following order:\n");
    printf("\t                                            1. Adjacency\n");
    printf("\t                                            2. Base Counting\n");
    printf("\t                                            3. Shouji\n");
    printf("\t                                            4. Q-Gram (short)\n");
    printf("\t                                            5. MAGNET\n");
    printf("\t                                            6. Hamming Distance\n");
    printf("\t                                            7. GateKeeper\n");
    printf("\t                                            8. SneakySnake\n");
    printf("\t                                            9. GRIM\n");
    printf("\t                                            10. KSW2\n");
    printf("\t                                            11. Edlib\n");
    printf("\t                                            12. GRIM\n");
    printf("\t                                            13. Q-Gram (long)\n");
    printf("\t--full-edlib                              Run Edlib without a limit in the edit distance threshold\n");
    printf("\t--help                                    Show this message\n");
}

int main(int argc, char **argv) {

    int active[FILTER_COUNT] = {0};
    int accepted[FILTER_COUNT] = {0};
    int false_positives[FILTER_COUNT] = {0};
    int false_negatives[FILTER_COUNT] = {0};
    double time_spent[FILTER_COUNT] = {0.0};
    double time[FILTER_COUNT] = {0.0};
    clock_t begin[FILTER_COUNT] = {0};
    clock_t end[FILTER_COUNT] = {0};

    char * filter_names[FILTER_COUNT] = {
            "Adjacency",
            "Base Counting",
            "Shouji",
            "Q-Gram (short)",
            "MAGNET",
            "Hamming Distance",
            "GateKeeper",
            "SneakySnake",
            "GRIM",
            "KSW2",
            "Edlib",
            "Grim (original)",
            "Q-Gram (long)"
    };

    for (int i = 0; i < FILTER_COUNT; i++) {
        active[i] = 1;
    }

    char * tss_filename;
    char * grid_filename = "grid";
    char * histogram_filename = "histogram";
    char * runtime_filename = "runtime.csv";

    int debug_mode = 0;
	int shouji_grid_size = 4;
	int adjacency_kmer_length = 5;
	int read_length = 1000; 
    int read_count = 100;
    int min_edit_distance = 0;
    int max_edit_distance = 10;

    int opt;
    int opt_index = 0;
    const char * opt_string = "d:s:k:m:f:c:i:g:h:r:u:l:a:t:p:";
    while((opt = getopt_long(argc, argv, opt_string, long_options, &opt_index)) >= 0)
    {
        switch(opt)
        {
            case 'd': // histogram filename
                debug_mode = atoi(optarg);
                break;
            case 's':
                shouji_grid_size = atoi(optarg);
                break;
            case 'k': 
                adjacency_kmer_length = atoi(optarg);
                break;
            case 'm': //set activation for each filter individually
                read_length = atoi(optarg);
                break;
            case 'f': 
                tss_filename = optarg;
                break;
            case 'c':
                read_count = atoi(optarg);
                break;
            case 'g':
                grid_filename = optarg;
                break;
            case 'h':
                histogram_filename = optarg;
                break;
            case 'r':
                runtime_filename = optarg;
                break;
            case 'u': //set high end for edit distance threshold
                max_edit_distance = atoi(optarg);
                break;
            case 'l': //set low end for edit distance threshold
                min_edit_distance = atoi(optarg);
                break;
            case 'a':
                if (strlen(optarg) == FILTER_COUNT) {
                    for (int i = 0; i < FILTER_COUNT; i++) {
                        active[i] = optarg[i] - '0';
                    }
                } else {
                    fprintf(stderr, "activation or deactivation values must be given for all filters:\n");
                    for (int i = 0; i < FILTER_COUNT-1; i++) {
                        fprintf(stderr, "\t%i. %s\n", i+1, filter_names[i]);
                    }
                    exit(1);
                }
                break;
            case 256:
                print_help(debug_mode, shouji_grid_size, adjacency_kmer_length, read_length, tss_filename, read_count, grid_filename, histogram_filename, runtime_filename, max_edit_distance, min_edit_distance);
                exit(0);
                break;
            case '?':
                fprintf(stderr, "unknown option: %c\n", optopt);
                print_help(debug_mode, shouji_grid_size, adjacency_kmer_length, read_length, tss_filename, read_count, grid_filename, histogram_filename, runtime_filename, max_edit_distance, min_edit_distance);
                exit(1);
                break;
        }
    }

    if (tss_filename == NULL) {
        fprintf(stderr, "-f is a required option. Please specify a .tss filename.\n");
        print_help(debug_mode, shouji_grid_size, adjacency_kmer_length, read_length, tss_filename, read_count, grid_filename, histogram_filename, runtime_filename, max_edit_distance, min_edit_distance);
        exit(1); 
    }

    /*
     * Read length is the MAX read length
     */

    int max_read_length = read_length;
	int n;
	FILE * tss_fp;
	char * line = NULL;
	size_t len = 0;
	ssize_t read;
	char *p;
	int j,i;
	int loop_edit_distance;
	int nw_accepted;
	int error_threshold;

    char * ref_seq = 0;
    ref_seq = (char *) calloc (max_read_length, sizeof(char));
    char * read_seq = 0;
    read_seq = (char *) calloc (max_read_length, sizeof(char));

    if (ref_seq == NULL) {
        fprintf(stderr, "ref_seq calloc failed\n");
        exit(1);
    }

    if (read_seq == NULL) {
        fprintf(stderr, "read_seq calloc failed\n");
        exit(1);
    }
    //generate csv header
    char header[512];
    csv_header(header, filter_names, 0);

    //csv accept/reject "grid" output file
    FILE * grid_fp;
    const char * directory = "results/";
    const char * filetype = ".csv";
    char grid_filename_buffer[512];

    mkdir(directory, 0777);

    // PRINT EXECUTION TIMES TO CSV
    FILE * runtime_fp;

    char runtime_filename_buffer[512];
    sprintf(runtime_filename_buffer, "%s%s", directory, runtime_filename);
    runtime_fp = fopen(runtime_filename_buffer, "w");

    fprintf(runtime_fp, "threshold,%s\n", header);

    //DATA STRUCTURE FOR HISTOGRAM
    int * histogram = (int *) calloc(100+1, sizeof(int));

    //loop_edit_distance gives edit distance threshold in percent... iterates over whole percentage points
	for (loop_edit_distance = min_edit_distance; loop_edit_distance<=max_edit_distance; loop_edit_distance++) {
		printf("<-------------------Levenshtein Distance = %d%%------------------->\n", loop_edit_distance);

        //print results to grid
        sprintf(grid_filename_buffer, "%s%s_%d%s", directory, grid_filename, loop_edit_distance, filetype);
        grid_fp = fopen(grid_filename_buffer, "w");
        fprintf(grid_fp, "%s,%s\n", header, GRID_CSV_SUFFIX); //TODO: modify to match correct output sequence

        memset(false_positives, 0, sizeof(false_positives));
        memset(false_negatives, 0, sizeof(false_negatives));
        memset(time, 0, sizeof(time));
        memset(begin, 0, sizeof(begin));
        memset(end, 0, sizeof(end));

        tss_fp = fopen(tss_filename, "r");

        char gridBuffer[GRID_LINE_LENGTH*GRID_LINES];
        char gridLineBuffer[GRID_LINE_LENGTH];
        char *gridp = gridBuffer; // maintains current position of the string in the buffer
        gridBuffer[0] = '\0';

        for (i = 1; i <= read_count; i++) {
            read = getline(&line, &len, tss_fp);
            j = 1;
            for (p = strtok(line, "\t"); p != NULL; p = strtok(NULL, "\t")) {
                // Read length is strlen(p)
                if (j == 1) {
                    read_length = strlen(p);
                    if (read_length > max_read_length) { // we truncate when above the max read length
                        read_length = max_read_length;
                    }
                    for (n = 0; n < read_length; n++) {
                        read_seq[n] = p[n];
                    }
                } else if (j == 2) {
                    if (strlen(p) < read_length) {
                        read_length = strlen(p); //TODO: check if this call is expensive
                    }
                    for (n = 0; n < read_length; n++) {
                        ref_seq[n] = p[n];
                    }
                }
                j = j + 1;
            }

            error_threshold = (loop_edit_distance * read_length) / 100;
            
            memset(accepted, 0, sizeof(accepted));

            // Adjacency Filter
            if (active[ADJACENCY_INDEX]) {
                begin[ADJACENCY_INDEX] = clock();
                accepted[ADJACENCY_INDEX] = AdjacencyFilter(read_length, ref_seq, read_seq, error_threshold, adjacency_kmer_length, debug_mode); //SlidingWindow(ReadLength, ref_seq, read_seq, error_threshold, shouji_grid_size, DebugMode);
                end[ADJACENCY_INDEX] = clock();
            }

            // Base Counting
            if (active[BASE_COUNTING_INDEX]) {
                begin[BASE_COUNTING_INDEX] = clock();
                accepted[BASE_COUNTING_INDEX] = baseCounting(read_length, ref_seq, read_seq, error_threshold, debug_mode);
                end[BASE_COUNTING_INDEX] = clock();
            }

            // Shouji
            if (active[SHOUJI_INDEX]) {
                begin[SHOUJI_INDEX] = clock();
                accepted[SHOUJI_INDEX] = Shouji(read_length, ref_seq, read_seq, error_threshold, shouji_grid_size, debug_mode);
                end[SHOUJI_INDEX] = clock();
            }

            // QGram
            if (active[QGRAM_INDEX]) {
                begin[QGRAM_INDEX] = clock();
                accepted[QGRAM_INDEX] = qgram(read_length, ref_seq, read_seq, error_threshold, 5);
                end[QGRAM_INDEX] = clock();
            }

            //MAGNET
            if (active[MAGNET_INDEX]) {
                begin[MAGNET_INDEX] = clock();
                accepted[MAGNET_INDEX] = MAGNET_DC(read_length, ref_seq, read_seq, error_threshold, debug_mode);
                end[MAGNET_INDEX] = clock();
            }

            //Hamming Distance
            if (active[HD_INDEX]) {
                begin[HD_INDEX] = clock();
                accepted[HD_INDEX] = HD(read_length, ref_seq, read_seq, error_threshold, debug_mode);
                end[HD_INDEX] = clock();
            }

            //SHD
            if (active[SHD_INDEX]) {
                begin[SHD_INDEX] = clock();
                accepted[SHD_INDEX] = SHD(read_length, ref_seq, read_seq, error_threshold, debug_mode);
                end[SHD_INDEX] = clock();
            }

            //SneakySnake
            if (active[SNEAKYSNAKE_INDEX]) {
                int SSWindow = read_length;
//                if (read_length > 5000) {
//                    SSWindow = error_threshold + 500;
//                }
                begin[SNEAKYSNAKE_INDEX] = clock();
                accepted[SNEAKYSNAKE_INDEX] = SneakySnake(read_length, ref_seq, read_seq, error_threshold, SSWindow, debug_mode, SSWindow);
                end[SNEAKYSNAKE_INDEX] = clock();
            }

            //GRIM Filter with updated threshold
            if (active[GRIM_INDEX]) {
                begin[GRIM_INDEX] = clock();
                accepted[GRIM_INDEX] = grim_original_tweak(read_length, ref_seq, read_seq, error_threshold, 5);
                end[GRIM_INDEX] = clock();
            }

            //KSW2
            if (active[KSW2_INDEX]) {
                begin[KSW2_INDEX] = clock();
                align(ref_seq, read_seq, 1, -2, 2, 1);
                end[KSW2_INDEX] = clock();
                accepted[KSW2_INDEX] = 1; // set KSW2 to accept everything
            }

            if (active[GRIM_ORIGINAL_INDEX]) {
                begin[GRIM_ORIGINAL_INDEX] = clock();
                accepted[GRIM_ORIGINAL_INDEX] = grim_original(read_length, ref_seq, read_seq, error_threshold, 5);
                end[GRIM_ORIGINAL_INDEX] = clock();
            }

            if (active[QGRAM_LONG_INDEX]) {
                begin[QGRAM_LONG_INDEX] = clock();
                accepted[QGRAM_LONG_INDEX] = qgram_hash(read_length, ref_seq, read_seq, error_threshold, 15);
                end[QGRAM_LONG_INDEX] = clock();
            }

            //Edlib
            EdlibAlignResult resultEdlib;
            if (full_edlib && loop_edit_distance == max_edit_distance) { // only do full edlib on the last run or if short reads
                begin[EDLIB_INDEX] = clock();
                resultEdlib = edlibAlign(ref_seq, read_length, read_seq, read_length, edlibNewAlignConfig(read_length, EDLIB_MODE_NW, EDLIB_TASK_DISTANCE, NULL, 0));
                end[EDLIB_INDEX] = clock();
            } else {
                begin[EDLIB_INDEX] = clock();
                resultEdlib = edlibAlign(ref_seq, read_length, read_seq, read_length, edlibNewAlignConfig(error_threshold, EDLIB_MODE_NW, EDLIB_TASK_DISTANCE, NULL, 0));
                end[EDLIB_INDEX] = clock();
            }

			if (resultEdlib.editDistance==-1) {
                if (loop_edit_distance == max_edit_distance) {
                    histogram[100]++;
                }
                accepted[EDLIB_INDEX] = 0;
            } else {
                if (loop_edit_distance == max_edit_distance && resultEdlib.editDistance >= 0 && resultEdlib.editDistance <= read_length) {
                    histogram[resultEdlib.editDistance*100/read_length]++;
                }
                accepted[EDLIB_INDEX] = (resultEdlib.editDistance <= error_threshold);
            }
			edlibFreeAlignResult(resultEdlib);

            nw_accepted = accepted[EDLIB_INDEX];

            //only print to file every 1024 lines (i starts at 1)
            //we assume each line is less than 512 chararcters

            // TODO: automate this to adapt to FILTER_COUNT
            int pos = 0;
            for(int i = 0; i < FILTER_COUNT; i++) {
                int consumed = sprintf(&gridLineBuffer[pos], "%i,", accepted[i]);
                pos += consumed;
            }
            sprintf(&gridLineBuffer[pos], "%i,%i,%i\n", resultEdlib.editDistance, error_threshold, read_length);
            gridp = linear_strcat(gridp, gridLineBuffer);

            if (i % GRID_LINES == 0 || i == read_count) {
                fprintf(grid_fp, "%s", gridBuffer);
                memset(gridBuffer, 0, sizeof(gridBuffer));
                gridp = gridBuffer; // reset pointer to beginning of grid
                gridBuffer[0] = '\0';
            }

            for (int i = 0; i < FILTER_COUNT; i++) {
                if (accepted[i] == 0 && nw_accepted == 1) false_negatives[i]++;
                else if (accepted[i] == 1 && nw_accepted == 0) false_positives[i]++;
                time[i] += end[i] - begin[i];
            }
		}

        for (int i = 0; i < FILTER_COUNT; i++) {
            time_spent[i] = (double) time[i] / CLOCKS_PER_SEC;

            // if a filter exceeds the total time limit, we deactivate it for the rest of the benchmark
            if (time_spent[i] > TIME_LIMIT_IN_SECONDS) {
                active[i] = 0;
            }
        }

        fprintf(runtime_fp, "%i,", loop_edit_distance);
        for (int i = 0; i < FILTER_COUNT; i++) {
            if (i != FILTER_COUNT-1) fprintf(runtime_fp, "%5.6f,", time_spent[i]);
            else  fprintf(runtime_fp, "%5.4f\n", time_spent[i]);
        }

		fclose(tss_fp);
        fclose(grid_fp);
    }

    // close file for runtime output
    fclose(runtime_fp);

    free(read_seq);
    free(ref_seq);

    FILE * histogram_fp;

    const char * histogram_directory = "results/";
    // if (argc > 9) {
    //     histogram_filename = argv[9];
    // }
    const char * histogram_filetype = ".csv";
    char histogramFileNameBuffer[512];

    sprintf(histogramFileNameBuffer, "%s%s%s", histogram_directory, histogram_filename, histogram_filetype);
    histogram_fp = fopen(histogramFileNameBuffer, "w");

    for (int i = 0; i <= 100; i++) {
        fprintf(histogram_fp, "%d, ", i);
    }
    fprintf(histogram_fp, "\n");
    for (int i = 0; i <= 100; i++) {
        fprintf(histogram_fp, "%d, ", histogram[i]);
    }

    fclose(histogram_fp);
    free(histogram);

	return 0;
}


