/*
 * MPI_Tile.cpp
 *
 *  Created on: 2015/11/04
 *      Author: stomo
 */

#include <iostream>
#include <cstdlib>
#include <mpi.h>
#include "MPI_Tile.hpp"

void TileSend( BMatrix *M, const int recver_rank, const int tag )
{
	int data_size = M->m() * M->n();
	MPI::COMM_WORLD.Send ( M->top(), data_size, MPI::DOUBLE, recver_rank, tag );
}

void TileRecv( BMatrix *M, const int sender_rank, const int tag )
{
	int data_size = M->m() * M->n();
	MPI::COMM_WORLD.Recv ( M->top(), data_size, MPI::DOUBLE, sender_rank, tag );
}

void TileBcast( BMatrix *M, const int sender_rank )
{
	int data_size = M->m() * M->n();
	MPI::COMM_WORLD.Bcast ( M->top(), data_size, MPI::DOUBLE, sender_rank );
}

MPI::Request TileISend( BMatrix *M, const int recver_rank, const int tag )
{
	int data_size = M->m() * M->n();
	return MPI::COMM_WORLD.Isend ( M->top(), data_size, MPI::DOUBLE, recver_rank, tag );
}

MPI::Request TileIRecv( BMatrix *M, const int sender_rank, const int tag )
{
	int data_size = M->m() * M->n();
	return MPI::COMM_WORLD.Irecv ( M->top(), data_size, MPI::DOUBLE, sender_rank, tag );
}

