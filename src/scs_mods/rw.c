#include "rw.h"
#include "scs.h"
#include "util.h"
#include "amatrix.h"
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>

/* writes/reads problem data to/from filename */
/* This is a VERY naive implementation, doesn't care about portability etc */

static void write_scs_cone(const ScsCone *k, FILE *fout) {
  fwrite(&(k->f), sizeof(scs_int), 1, fout);
  fwrite(&(k->l), sizeof(scs_int), 1, fout);
  fwrite(&(k->qsize), sizeof(scs_int), 1, fout);
  fwrite(k->q, sizeof(scs_int), k->qsize, fout);
  fwrite(&(k->ssize), sizeof(scs_int), 1, fout);
  fwrite(k->s, sizeof(scs_int), k->ssize, fout);
  fwrite(&(k->ep), sizeof(scs_int), 1, fout);
  fwrite(&(k->ed), sizeof(scs_int), 1, fout);
  fwrite(&(k->psize), sizeof(scs_int), 1, fout);
  fwrite(k->p, sizeof(scs_float), k->psize, fout);
}

static ScsCone * read_scs_cone(FILE *fin) {
  size_t bytes;
  ScsCone * k = (ScsCone *)scs_calloc(1, sizeof(ScsCone));
  bytes = fread(&(k->f), sizeof(scs_int), 1, fin);
  if (!bytes) scs_printf("no data read!\n");
  bytes = fread(&(k->l), sizeof(scs_int), 1, fin);
  if (!bytes) scs_printf("no data read!\n");
  bytes = fread(&(k->qsize), sizeof(scs_int), 1, fin);
  if (!bytes) scs_printf("no data read!\n");
  k->q = scs_calloc(k->qsize, sizeof(scs_int));
  bytes = fread(k->q, sizeof(scs_int), k->qsize, fin);
  if (!bytes) scs_printf("no data read!\n");
  bytes = fread(&(k->ssize), sizeof(scs_int), 1, fin);
  if (!bytes) scs_printf("no data read!\n");
  k->s = scs_calloc(k->ssize, sizeof(scs_int));
  bytes = fread(k->s, sizeof(scs_int), k->ssize, fin);
  if (!bytes) scs_printf("no data read!\n");
  bytes = fread(&(k->ep), sizeof(scs_int), 1, fin);
  if (!bytes) scs_printf("no data read!\n");
  bytes = fread(&(k->ed), sizeof(scs_int), 1, fin);
  if (!bytes) scs_printf("no data read!\n");
  bytes = fread(&(k->psize), sizeof(scs_int), 1, fin);
  if (!bytes) scs_printf("no data read!\n");
  k->p = scs_calloc(k->psize, sizeof(scs_int));
  bytes = fread(k->p, sizeof(scs_float), k->psize, fin);
  if (!bytes) scs_printf("no data read!\n");
  RETURN k;
}

static void write_scs_stgs(const ScsSettings *s, FILE *fout) {
  /* Warm start to false for now */
  scs_int warm_start = 0;
  fwrite(&(s->normalize), sizeof(scs_int), 1, fout);
  fwrite(&(s->scale), sizeof(scs_float), 1, fout);
  fwrite(&(s->rho_x), sizeof(scs_float), 1, fout);
  fwrite(&(s->max_iters), sizeof(scs_int), 1, fout);
  fwrite(&(s->eps), sizeof(scs_float), 1, fout);
  fwrite(&(s->alpha), sizeof(scs_float), 1, fout);
  fwrite(&(s->cg_rate), sizeof(scs_float), 1, fout);
  fwrite(&(s->verbose), sizeof(scs_int), 1, fout);
  fwrite(&warm_start, sizeof(scs_int), 1, fout);
  fwrite(&(s->acceleration_lookback), sizeof(scs_int), 1, fout);
  /* Do not write the write_data_filename */
}

static ScsSettings * read_scs_stgs(FILE *fin) {
  size_t bytes;
  ScsSettings * s = (ScsSettings *)scs_calloc(1, sizeof(ScsSettings));
  bytes = fread(&(s->normalize), sizeof(scs_int), 1, fin);
  if (!bytes) scs_printf("no data read!\n");
  bytes = fread(&(s->scale), sizeof(scs_float), 1, fin);
  if (!bytes) scs_printf("no data read!\n");
  bytes = fread(&(s->rho_x), sizeof(scs_float), 1, fin);
  if (!bytes) scs_printf("no data read!\n");
  bytes = fread(&(s->max_iters), sizeof(scs_int), 1, fin);
  if (!bytes) scs_printf("no data read!\n");
  bytes = fread(&(s->eps), sizeof(scs_float), 1, fin);
  if (!bytes) scs_printf("no data read!\n");
  bytes = fread(&(s->alpha), sizeof(scs_float), 1, fin);
  if (!bytes) scs_printf("no data read!\n");
  bytes = fread(&(s->cg_rate), sizeof(scs_float), 1, fin);
  if (!bytes) scs_printf("no data read!\n");
  bytes = fread(&(s->verbose), sizeof(scs_int), 1, fin);
  if (!bytes) scs_printf("no data read!\n");
  bytes = fread(&(s->warm_start), sizeof(scs_int), 1, fin);
  if (!bytes) scs_printf("no data read!\n");
  bytes = fread(&(s->acceleration_lookback), sizeof(scs_int), 1, fin);
  if (!bytes) scs_printf("no data read!\n");
  RETURN s;
}

static void write_amatrix(const ScsMatrix * A, FILE * fout) {
  scs_int Anz = A->p[A->n];
  fwrite(&(A->m), sizeof(scs_int), 1, fout);
  fwrite(&(A->n), sizeof(scs_int), 1, fout);
  fwrite(A->p, sizeof(scs_int), A->n + 1, fout);
  fwrite(A->x, sizeof(scs_float), Anz, fout);
  fwrite(A->i, sizeof(scs_int), Anz, fout);
}

static ScsMatrix * read_amatrix(FILE * fin) {
  scs_int Anz;
  size_t bytes;
  ScsMatrix * A = (ScsMatrix *)scs_calloc(1, sizeof(ScsMatrix));
  bytes = fread(&(A->m), sizeof(scs_int), 1, fin);
  if (!bytes) scs_printf("no data read!\n");
  bytes = fread(&(A->n), sizeof(scs_int), 1, fin);
  if (!bytes) scs_printf("no data read!\n");
  A->p = scs_calloc(A->n + 1, sizeof(scs_int));
  bytes = fread(A->p, sizeof(scs_int), A->n + 1, fin);
  if (!bytes) scs_printf("no data read!\n");
  Anz = A->p[A->n];
  A->x = scs_calloc(Anz, sizeof(scs_float));
  A->i = scs_calloc(Anz, sizeof(scs_int));
  bytes = fread(A->x, sizeof(scs_float), Anz, fin);
  if (!bytes) scs_printf("no data read!\n");
  bytes = fread(A->i, sizeof(scs_int), Anz, fin);
  if (!bytes) scs_printf("no data read!\n");
  RETURN A;
}

static void write_scs_data(const ScsData *d, FILE *fout) {
  fwrite(&(d->m), sizeof(scs_int), 1, fout);
  fwrite(&(d->n), sizeof(scs_int), 1, fout);
  fwrite(d->b, sizeof(scs_float), d->m, fout);
  fwrite(d->c, sizeof(scs_float), d->n, fout);
  write_scs_stgs(d->stgs, fout);
  write_amatrix(d->A, fout);
}

static ScsData * read_scs_data(FILE *fin) {
  size_t bytes;
  ScsData * d = (ScsData *)scs_calloc(1, sizeof(ScsData));
  bytes = fread(&(d->m), sizeof(scs_int), 1, fin);
  if (!bytes) scs_printf("no data read!\n");
  bytes = fread(&(d->n), sizeof(scs_int), 1, fin);
  if (!bytes) scs_printf("no data read!\n");
  d->b = scs_calloc(d->m, sizeof(scs_float));
  d->c = scs_calloc(d->n, sizeof(scs_float));
  bytes = fread(d->b, sizeof(scs_float), d->m, fin);
  if (!bytes) scs_printf("no data read!\n");
  bytes = fread(d->c, sizeof(scs_float), d->n, fin);
  if (!bytes) scs_printf("no data read!\n");
  d->stgs = read_scs_stgs(fin);
  d->A = read_amatrix(fin);
  RETURN d;
}

void SCS(write_data)(const ScsData *d, const ScsCone *k) {
  FILE* fout = fopen(d->stgs->write_data_filename, "wb");
  uint32_t scs_int_sz = (uint32_t)sizeof(scs_int);
  uint32_t scs_float_sz = (uint32_t)sizeof(scs_float);
  scs_printf("writing data to %s\n", d->stgs->write_data_filename);
  fwrite(&(scs_int_sz), sizeof(uint32_t), 1, fout);
  fwrite(&(scs_float_sz), sizeof(uint32_t), 1, fout);
  write_scs_cone(k, fout);
  write_scs_data(d, fout);
  fclose(fout);
}

scs_int SCS(read_data)(const char * filename, ScsData ** d, ScsCone ** k) {
  uint32_t file_int_sz;
  uint32_t file_float_sz;
  size_t bytes;
  FILE* fin = fopen(filename, "rb");
  if (!fin) {
    scs_printf("Error reading file %s\n", filename);
    RETURN -1;
  }
  scs_printf("Reading data from %s\n", filename);
  bytes = fread(&(file_int_sz), sizeof(uint32_t), 1, fin);
  if (!bytes) scs_printf("no data read!\n");
  bytes = fread(&(file_float_sz), sizeof(uint32_t), 1, fin);
  if (!bytes) scs_printf("no data read!\n");
  
  if (file_int_sz != (uint32_t)sizeof(scs_int)) {
    scs_printf("Error, sizeof(file int) is %lu, but scs expects sizeof(int) "
        "%lu, scs should be recompiled with correct flags.\n",
        (unsigned long)file_int_sz, (unsigned long)sizeof(scs_int));
    fclose(fin);
    RETURN -1;
  }
  if (file_float_sz != (uint32_t)sizeof(scs_float )) {
    scs_printf("Error, sizeof(file float) is %lu, but scs expects "
        "sizeof(float) %lu, scs should be recompiled with the correct flags.\n",
        (unsigned long)file_float_sz, (unsigned long)sizeof(scs_float));
    fclose(fin);
    RETURN -1;
  }

  *k = read_scs_cone(fin);
  *d = read_scs_data(fin);
  fclose(fin);
  RETURN 0;
}
