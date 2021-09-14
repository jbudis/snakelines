/* Philip T.L.C. Clausen Jan 2017 plan@dtu.dk */

/*
 * Copyright (c) 2017, Philip Clausen, Technical University of Denmark
 * All rights reserved.
 * 
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 * 
 *		http://www.apache.org/licenses/LICENSE-2.0
 * 
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
*/
#define _XOPEN_SOURCE 600
//cc -Wall -O3 -std=c99 -c -o seqmenttree seqmenttree.c
#ifndef SEQMENTTREE
#define SEQMENTTREE 1
typedef struct seqmentTree SeqmentTree;
typedef struct seqmentTrees SeqmentTrees;
struct seqmentTree {
	unsigned n;
	unsigned size;
	struct seqmentTrees *root;
};
struct seqmentTrees {
	unsigned start;
	unsigned end;
	unsigned covered;
	struct seqmentTrees *branch[2];
};
#endif

SeqmentTree * initializeSeqmentTree(long unsigned size);
SeqmentTree * initSeqmentTree(SeqmentTree *src, const unsigned start, const unsigned end);
void reallocSeqmentTree(SeqmentTree *src);
unsigned addSeqmentTrees(SeqmentTrees *root, SeqmentTrees *node);
int growSeqmentTree(SeqmentTree *src, const unsigned start, const unsigned end);
unsigned queSeqmentTree(SeqmentTrees *src, const unsigned start, const unsigned end);
