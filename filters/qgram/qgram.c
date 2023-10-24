//
// Created by mdr on 06/07/2021.
//

#include "qgram.h"
#include <stdlib.h>
#include <string.h>
#include "../../libs/khash.h"
#include <stdio.h>

KHASH_MAP_INIT_INT(m32, char)

static int index_from_qgram(char * gram, int qGramLength) { // calculate natural order index from QGram
    unsigned int index = 0;
    unsigned int offset = 1 << 2 * (qGramLength - 1); // is equal to 4^qGramLength / 4
    for (int i = 0; i < qGramLength; i++) {
        if (gram[i] == 'T') {
            index += offset;
        } else if (gram[i] == 'G') {
            index += 2 * offset;
        } else if (gram[i] == 'C') {
            index += 3 * offset;
        } // we do not offset the index for A... this saves a case
        offset /= 4;
    }
    return index;
}

static int getOccurrenceTableLength(int qGramLength) {
    return 1 << 2*qGramLength;
}


// FLASH
int qgram(int ReadLength, const char RefSeq[], const char ReadSeq[], int ErrorThreshold, int qGramLength) {
    // assumes ReadLength == RefLength

    int occurrenceTableLength = getOccurrenceTableLength(qGramLength);
    int * occurrenceCount = (int *) calloc(occurrenceTableLength, sizeof(int));

    char refGram[qGramLength];
    char readGram[qGramLength];
    int qGramCount = ReadLength - qGramLength + 1;

    for (int i = 0; i < qGramCount; i++) {
        memcpy(refGram, &RefSeq[i], qGramLength);
        occurrenceCount[index_from_qgram(refGram, qGramLength)]++;
    }

    int threshold = qGramLength*ErrorThreshold;

    int count = 0;
    for (int i = 0; i < qGramCount; i++) {
        memcpy(readGram, &ReadSeq[i], qGramLength);
        if (occurrenceCount[index_from_qgram(readGram, qGramLength)] > 0) {
            occurrenceCount[index_from_qgram(readGram, qGramLength)]--;
        } else {
            count++;
            if (count > threshold) {
                break;
            }
        }
    }

    free(occurrenceCount);

    return count <= threshold;
}

void update_hash_table(khash_t(m32) * h, int idx, int val) {
    int ret, is_missing;
    khint_t k = kh_get(m32, h, idx);          // query the hash table
    is_missing = (k == kh_end(h));   // test if the key is present
    if (is_missing) {
        k = kh_put(m32, h, idx, &ret);     // insert a key to the hash table
        if (!ret) kh_del(m32, h, k);
        kh_value(h, k) = val;             // set the value
    } else {
        kh_val(h, k) += val;
    }
}

int qgram_hash(int ReadLength, const char RefSeq[], const char ReadSeq[], int ErrorThreshold, int qGramLength) {
    // assumes ReadLength == RefLength
    int is_missing;

    khash_t(m32) *h = kh_init(m32);  // allocate a hash table

    char refGram[qGramLength];
    char readGram[qGramLength];
    int qGramCount = ReadLength - qGramLength + 1;
    int threshold = qGramLength*ErrorThreshold;
    int idx = 0;

    for (int i = 0; i < qGramCount; i++) {
        memcpy(refGram, &RefSeq[i], qGramLength);
        idx = index_from_qgram(refGram, qGramLength);
        update_hash_table(h, idx, 1); //increment the hash table
    }

    khint_t k;

    int count = 0;
    for (int i = 0; i < qGramCount; i++) {
        memcpy(readGram, &ReadSeq[i], qGramLength);
        idx = index_from_qgram(readGram, qGramLength);

        k = kh_get(m32, h, idx);          // query the hash table
        is_missing = (k == kh_end(h));   // test if the key is present

        if (is_missing) {
            count++;
            if (count > threshold) {
                break;
            }
        } else {
            update_hash_table(h, idx, -1);
        }
    }

    // traverse table to complete sum
    for (k = kh_begin(h); k != kh_end(h); ++k)  // traverse
        if (kh_exist(h, k))          // test if a bucket contains data
            count += abs(kh_val(h, k));

    kh_destroy_m32(h);
    return count <= threshold;
}


//
//int qgram_hash(int ReadLength, const char RefSeq[], const char ReadSeq[], int ErrorThreshold, int qGramLength) {
//    int qGramCount = ReadLength - qGramLength;
//    khash_t(occ_hash) *occ_table = kh_init(occ_hash);
//
//    char refGram[qGramLength];
//    char readGram[qGramLength];
//
//    int absent, idx;
//
//    for (int i = 0; i < qGramCount; i++) {
//
//        memcpy(readGram, &ReadSeq[i], qGramLength);
//
//        idx = index_from_qgram(readGram, qGramLength);
//
//        khint_t k;
//
//        k = kh_get(occ_hash, occ_table, idx);
//
//        if (k == kh_end(occ_table)) { // key is not present in hash table
//            khint_t entry = kh_put(occ_hash, occ_table, idx, &absent);
////            printf("absent: %i", absent);
//            kh_val(occ_table, entry) = 1;
//        } else {
//            kh_val(occ_table, k) += 1;
//        }
//    }
//
//    for (int i = 0; i < qGramCount; i++) {
//
//        memcpy(refGram, &RefSeq[i], qGramLength);
//
//        idx = indexFromQGram(refGram, qGramLength);
//
//        khint_t k;
//
//        k = kh_get(occ_hash, occ_table, idx);
//
//        if (k == kh_end(occ_table)) { // key is not present in hash table
//            khint_t entry = kh_put(occ_hash, occ_table, idx, &absent);
//            kh_val(occ_table, entry) = -1;
//        } else {
//            kh_val(occ_table, k) -= 1;
//        }
//    }
//
//    int errorTotal = 0;
//    for (khint_t k = kh_begin(occ_table); k != kh_end(occ_table); ++k) {
//       errorTotal += abs(kh_val(occ_table, k));
//    }
//
//    kh_destroy(occ_hash, occ_table);
//
//    return errorTotal <= 2*qGramLength*ErrorThreshold;
//}

int qgram_avx2(int ReadLength, const char RefSeq[], const char ReadSeq[], int ErrorThreshold, int qGramLength) {
    return 0;
}