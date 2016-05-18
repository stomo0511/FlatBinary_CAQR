//============================================================================
// Name        : MPI_LL.cpp
// Author      : T. Suzuki
// Version     : isend/irecv
// Copyright   :
// Description :
//  Date       : 2015/11/09
//============================================================================

#include <cstdlib>
#include <cassert>
#include <mpi.h>
#include <omp.h>

#include <plasma.h>
#include <cblas.h>
#include <lapacke.h>

#include <Matrix.hpp>
#include <TMatrix.hpp>
#include <CoreBlasTile.hpp>
#include <Progress.hpp>
#include "MPI_Tile.hpp"

#ifdef VTRACE
#include <vt_user.h>
#endif

using namespace std;

int main(int argc, char * argv[])
{
	#ifdef VTRACE
	VT_OFF();
	#endif

	if (argc < 5)
    {
		cerr << "Usage: LL [M] [N] [NB] [IB]\n";
		exit(EXIT_FAILURE);
    }

	const int M  = atoi(argv[1]); // n. of rows of the matrix
	const int N  = atoi(argv[2]); // n. of columns of the matrix
	const int NB = atoi(argv[3]); // tile size
	const int IB = atoi(argv[4]); // inner blocking size

	assert( M >= N );
	assert( NB % IB == 0 );

	// MPI initialize
	MPI::Init ();
//	int req = MPI_THREAD_SERIALIZED;
//	int prv = MPI::Init_thread(req);
//	cout << "req = " << req << ", PRV = " << prv << endl;

	int my_rank = MPI::COMM_WORLD.Get_rank ( );
	int n_procs = MPI::COMM_WORLD.Get_size ( );

//	if (n_procs == 1)
//	{
//		cerr << "There is only one MPI process\n";
//		MPI::Finalize();
//		exit(EXIT_FAILURE);
//	}

	#ifdef DEBUG
	if (my_rank == 0)
    {
		cout << "M = " << M << ", N = " << N;
		cout << ", NB = " << NB << ", IB = " << IB << endl;
    }
	#endif

    ////////////////////////////////////////////////////////////////////////////
	// Matrix "A"
	TMatrix A( M, N, NB, NB, IB );

	// Initialize matrix "A"
	A.Set_Rnd( 20151030 );

	const int MT = A.mt();  // # tile rows
	const int NT = A.nt();  // # tile columns

	if (NT < n_procs)
    {
		cerr << "Too much processes\n";
		MPI_Finalize();
		return EXIT_FAILURE;
    }

	// Matrix "T"
	TMatrix T( MT*IB, NT*NB, IB, NB, IB );

	// Progress Table "Pt"
	Progress_Table Pt( MT, NT, NT );

    ////////////////////////////////////////////////////////////////////////////

	#ifdef DEBUG
	// Copy "A" to "cA" for debug
	Matrix cA( M, N );

	if (my_rank == 0)
    {
		A.Mat_Copy(cA);
    }
	#endif

	MPI::COMM_WORLD.Barrier();

	// Timer start
	double time;
	if (my_rank ==0)
		time = MPI::Wtime();

	#ifdef VTRACE
	VT_ON();
	#endif

	omp_set_num_threads(4);

//	#ifdef DEBUG
//	cout << "# of threads=" << omp_get_max_threads() << endl;
//	#endif

	#pragma omp parallel for schedule(dynamic,1)
	for ( int k=my_rank - n_procs; k<NT; k+=n_procs )
	{
		// Communication thread
		if (k == my_rank - n_procs)
		{
			double* Buf = new double [ (NB+IB)*NB ];

			for ( int ck=0; ck<NT-1; ck++ )
			{
				for ( int ci=ck; ci<MT; ci++ )
				{
					int sender = ck % n_procs;

					/////////////////////////////////
					// I'm a Sender
					if ( sender == my_rank )
					{
						Pt.check_waitIJK( ci,ck,ck );

						for ( int p = 0; p < n_procs; p++ )
						{
							int recver = p % n_procs;

							if ( recver != my_rank )
							{
								cblas_dcopy( NB*NB, A(ci,ck)->top(), 1, Buf, 1 );
								cblas_dcopy( IB*NB, T(ci,ck)->top(), 1, Buf+(NB*NB), 1 );

								MPI::COMM_WORLD.Send ( Buf, (NB+IB)*NB, MPI::DOUBLE, recver, TAG_G_A );
							}
						}
					}
					/////////////////////////////////
					// I'm a Receiver
					else
					{
						MPI::COMM_WORLD.Recv ( Buf, (NB+IB)*NB, MPI::DOUBLE, sender, TAG_G_A );

						cblas_dcopy( NB*NB, Buf, 1, A(ci,ck)->top(), 1 );
						cblas_dcopy( IB*NB, Buf+(NB*NB), 1, T(ci,ck)->top(), 1 );

						Pt.setIJK( ci,ck,ck,DONE );
					}
					/////////////////////////////////
				} // ci-loop end
			} // ck-loop end

			delete [] Buf;
		}
		else
		{
			// Computation threads
			for ( int l=0; l<k; l++ )
			{
				////////////////////////////////////////////////////////////////////
				// LARFB
				Pt.check_waitIJK( l,l,l );
				{
					#ifdef VTRACE
					VT_TRACER("LARFB");
					#endif
					LARFB( PlasmaLeft, PlasmaTrans, A(l,l), T(l,l), A(l,k) );
				}
				#ifdef DEBUG
				cout << "LARFB (" << l << "," << k << "," << l << ") : ";
				cout << omp_get_thread_num() << " of " << my_rank << endl;
				#endif
				// LARFB end
				////////////////////////////////////////////////////////////////////

				////////////////////////////////////////////////////////////////////
				// SSRFB
				for ( int i=l+1; i<MT; i++ )
				{
					Pt.check_waitIJK( i,l,l );
					{
						#ifdef VTRACE
						VT_TRACER("SSRFB");
						#endif
						SSRFB( PlasmaLeft, PlasmaTrans, A(i,l), T(i,l), A(l,k), A(i,k) );
					}
					#ifdef DEBUG
					cout << "SSRFB (" << i << "," << k << "," << l << ") : ";
					cout << omp_get_thread_num() << " of " << my_rank << endl;
					#endif
				} // i-loop end
				// SSRFB end
				////////////////////////////////////////////////////////////////////
			} // l-loop end

			////////////////////////////////////////////////////////////////////
			// GEQRT
			{
				#ifdef VTRACE
					VT_TRACER("GEQRT");
					#endif
				GEQRT( A(k,k), T(k,k) );
			}
			Pt.setIJK(k,k,k,DONE); // Update progress table
			#ifdef DEBUG
			cout << "GEQRT (" << k << "," << k << "," << k << ") : ";
			cout << omp_get_thread_num() << " of " << my_rank << endl;
			#endif
			// GEQRT end
			////////////////////////////////////////////////////////////////////

			////////////////////////////////////////////////////////////////////
			// TSQRT
			for (int i=k+1; i<MT; i++)
			{
				{
					#ifdef VTRACE
						VT_TRACER("TSQRT");
						#endif
					TSQRT( A(k,k), A(i,k), T(i,k) );
				}
				Pt.setIJK(i,k,k,DONE);  // Update progress table
				#ifdef DEBUG
				cout << "TSQRT (" << i << "," << k << "," << k << ") : ";
				cout << omp_get_thread_num() << " of " << my_rank << endl;
				#endif
			} // i-loop end
			// TSQRT end
			////////////////////////////////////////////////////////////////////
		}
	}
	// Main loop (k-loop) end
	// omp parallel end

    // Timer stop
    MPI::COMM_WORLD.Barrier();
    if (my_rank ==0) {
    	time = MPI::Wtime() - time;
    	cout << M << ", " << time << endl;
    }

	#ifdef VTRACE
   VT_OFF();
   #endif

	#ifdef DEBUG
   for (int j=0; j<NT; j++)
    {
	   int iroot = j % n_procs;
	   for (int i=0; i<MT; i++)
		   TileBcast( A(i,j), iroot );
    	MPI::COMM_WORLD.Barrier();
     }
    for (int j=0; j<NT; j++)
     {
    	int iroot = j % n_procs;
    	for (int i=0; i<MT; i++)
    		TileBcast( T(i,j), iroot );
    	MPI::COMM_WORLD.Barrier();
     }

    if (my_rank == 0)
    {
        // Regenerate "Q"
        TMatrix Q(M,M,NB,NB,IB);

        // Set to the identity matrix
        Q.Set_Iden();

        // Make Orthogonal matrix "Q"
        dorgqr( A, T, Q );

        // Copy "Q" to "mQ"
        Matrix mQ( M, M );
        Q.Mat_Copy(mQ);
        // Regenerate "Q" End
        ////////////////////////////////////////////////////////////////////////

        ////////////////////////////////////////////////////////////////////////
        // Construct "R"
        Matrix mR( M, N );
        A.Mat_Copy(mR);

        // "R" is the upper triangler part of "A"
        for (int ii=0; ii<M; ii++)
        {
            for (int jj = 0; jj < N; jj++)
            {
                if ( ii > jj )
                    mR.Set_Val( ii, jj, (double)(0.0) );

            }
        }
        // Construct "R"
        ////////////////////////////////////////////////////////////////////////

        ////////////////////////////////////////////////////////////////////////
        // Check Orthogonarity
        // Set Id to the identity matrix
        Matrix Id( N, N );
        Id.Set_Iden();

        double alpha =  1.0;
        double beta  = -1.0;

        cblas_dsyrk(CblasColMajor, CblasUpper, CblasTrans,
                    N, M, alpha, mQ.top(), M, beta, Id.top(), N);

        char NRM = 'I';
        char UPL = 'U';

        double normQ = LAPACKE_dlansy( LAPACK_COL_MAJOR, NRM, UPL, N, Id.top(), N );
        cout << "norm(I-Q*Q') = " << normQ << endl;

        // Check Orthogonarity END
        ////////////////////////////////////////////////////////////////////////

        ////////////////////////////////////////////////////////////////////////
        // Check Residure Norm

        cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans,
                    M, N, M, alpha, mQ.top(), M, mR.top(), M, beta, cA.top(), M);

        normQ = LAPACKE_dlange( LAPACK_COL_MAJOR, NRM, M, N, cA.top(), M );
        cout << "norm(A-Q*R) = " << normQ << endl;

        // Check Residure END
        ///////////////////////////////////////////////////////////////////////
    }
    #endif

   MPI::Finalize();

   return EXIT_SUCCESS;
}
