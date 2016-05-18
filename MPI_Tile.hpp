/*
 * MPI_Tile.h
 *
 *  Created on: 2015/11/04
 *      Author: stomo
 */

#ifndef MPI_TILE_H_
#define MPI_TILE_H_

#define TAG_G_A 100
#define TAG_G_T 101
#define TAG_T_A 200
#define TAG_T_T 201
#define TAG_S_S 300

#include "BMatrix.hpp"
#include <mpi.h>

// Type declaration of the communication buffer for Transformation matrices
// MPI::Datatype t_mats;

void TileSend( const BMatrix *M, const int recver_rank, const int tag);
void TileRecv( BMatrix *M, const int sender_rank, const int tag);
void TileBcast( BMatrix *M, const int sender_rank);
MPI::Request TileISend( const BMatrix *M, const int recver_rank, const int tag );
MPI::Request TileIRecv( BMatrix *M, const int sender_rank, const int tag );

#endif /* MPI_TILE_H_ */
