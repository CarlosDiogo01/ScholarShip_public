/**************************************************************************
 *                                                                         *
 *         Java Grande Forum Benchmark Suite - Thread Version 1.0          *
 *                                                                         *
 *                            produced by                                  *
 *                                                                         *
 *                  Java Grande Benchmarking Project                       *
 *                                                                         *
 *                                at                                       *
 *                                                                         *
 *                Edinburgh Parallel Computing Centre                      *
 *                                                                         *
 *                email: epcc-javagrande@epcc.ed.ac.uk                     *
 *                                                                         *
 *      adapted from SciMark 2.0, author Roldan Pozo (pozo@cam.nist.gov)   *
 *                                                                         *
 *      This version copyright (c) The University of Edinburgh, 2001.      *
 *                         All rights reserved.                            *
 *                                                                         *
 **************************************************************************/

package sm.improved;

import java.util.Random;
import java.util.concurrent.BrokenBarrierException;
import java.util.concurrent.CyclicBarrier;
import java.util.concurrent.TimeUnit;
import java.lang.Object;

public class SOR
{

	public static double Gtotal = 0.0;
	public static final int cachelinesize = 128;
	public static volatile long sync[][];

	
	SOR(int total_threads){
		sync = new long  [total_threads][cachelinesize];
	}

	public final void SORrun(double omega, double G[][], int M, int N, int num_iterations, Random R) {		
		
		final int Mm1 			= M-1;
		final int Nm1 			= N-1;
			
		/* Creating Threads and Barrier */
		final int total_threads = JGFSORBench.nthreads;
		final CyclicBarrier barrier = new CyclicBarrier(total_threads);
		final Thread th[] 			= new Thread[total_threads-1];

		/* Measuring time with Guava StopWatch */
		StopWatch timerThreads	= StopWatch.createUnstarted();
		StopWatch timerMaster	= StopWatch.createUnstarted();
		
		final int iterations = (Mm1+2)/2;
		for (int i = 0; i < total_threads-1; i++){
			final int id = i+1;			
			th[i] = new Thread(new Runnable() {				
				public final void run(){
					
					int iter_per_thr	= iterations / total_threads;
					int iter_extra		= iterations % total_threads;
					
					if (id < iter_extra) {		//Threads that will receive the extra iterations (if any)
						iter_extra = 0;
						iter_per_thr++;
						
					}
					final int ilow		= (iter_per_thr * id + iter_extra) * 2+1; //begin = 1
					final int iupper	= ilow + iter_per_thr * 2;
					
					/* Applying First Touch - Using threads to Alloc and Initialize G */
					firstTouchAllocation(G, N, R, ilow-1, iupper-1);
					callBarrier(barrier);
					
					/* Start Measuring Time Working Threads */
					try {
						synchronized(timerThreads){
							timerThreads.start();
						}
					} catch (IllegalStateException e) {} // ignore start counter if watch was already started
					SORRunner(G, id, omega, num_iterations, ilow, iupper, Mm1, Nm1, total_threads);
				}
			});
			th[i].start();
		}
		
		int iter_per_thr = iterations / total_threads;
		
		//Threads that will receive the extra iterations (if any)
		if (0 < (iterations % total_threads)){
			iter_per_thr++;
		}
		final int iupper = iter_per_thr * 2;
		timerThreads.stop();
		
		/* Applying First Touch to MASTER */
		firstTouchAllocation(G, N, R, 0, iupper);
		callBarrier(barrier);

		/* Start Measuring Time MASTER */
		try {
			synchronized(timerMaster){
				timerMaster.start();
			}
		} catch (IllegalStateException e) {} // ignore start counter if watch was already started
		SORRunner(G, 0, omega, num_iterations, 1, iupper, Mm1, Nm1, total_threads);		
		killThreads(total_threads, th);
		timerMaster.stop();
		/******* End measuring time master *******/
		
		/* Computing time result */
		final long timeElapsedThreads	= timerThreads.elapsed(TimeUnit.MILLISECONDS);
		final long timeElapsedMaster	= timerMaster.elapsed(TimeUnit.MILLISECONDS);
		System.out.println(timeElapsedThreads + timeElapsedMaster);
		
		sumTotalG(G, Nm1);
	}	
	
	
	private void sumTotalG(double[][]G, final int Nm1){
		
		for (int i=1; i<Nm1; i++){
			for (int j=1; j<Nm1; j++){
				Gtotal += G[i][j];
			}
		}
	}
	
	private void killThreads(final int total_threads, final Thread[] th){
		
		for(int i=0; i<total_threads-1; i++){
			
			try 
			{
				th[i].join();
			}
			catch (InterruptedException e) {}
		}
	}
	
	private void callBarrier(final CyclicBarrier barrier){
		
		try
		{
			barrier.await();
		}
		catch (InterruptedException | BrokenBarrierException e)
		{
			e.printStackTrace();
		}
		
	}
	
	private void firstTouchAllocation(double[][] G, int N, Random R, final int ilow, final int iupper){
		
		for (int i=ilow; i < iupper; i++)
		{	
			/* Allocation in Parallel */
			G[i] = new double[N];
				
			/* Initialize in Parallel */
			for (int j=0; j<N; j++) {
				G[i][j] = R.nextDouble() * 1e-6;
			}
		}
	}
	
	
	
	public final void SORRunner(	final double G[][], final int id, final double omega, 
								final int num_iterations, final int ilow, final int iupper,
								final int Mm1, final int Nm1, final int total_threads)
	{
		final double omega_over_four 	= omega * 0.25;
		final double one_minus_omega 	= 1.0 - omega;
		final int num_iterations_2 		= 2 * num_iterations;
	
		// adjusting to remove the last iteration from the loop
		final int end = iupper - 2; 
		
		for (int p=0; p< num_iterations_2; p++) 
		{
			int begin = ilow + (p%2);
			
			// Dealing with the first row
			begin = first_row(begin, G, omega_over_four, one_minus_omega, Nm1); 
			
			// Dealing with the middle rows
			int i;
			for (i = begin; i < end; i += 2) 
			{
				middle_rows(G, omega_over_four, one_minus_omega, Nm1, i, G[i], G[i-1]);
			}
			
			// Dealing with the last iteration
			if (i == Mm1) 
			{
				last_row(G, omega_over_four, one_minus_omega, Nm1, i, G[i], G[i-1]);
			} 
			else if(i < Mm1)
			{
				middle_rows(G, omega_over_four, one_minus_omega, Nm1, i, G[i], G[i-1]);
			}

			call_barrier(id, total_threads);
		}

	}

	public int first_row(	int begin, final double[][] G, final double omega_over_four, 
							final double one_minus_omega, final int Nm1) {
		if(begin == 1) 
		{ 
			first_row(G, omega_over_four, one_minus_omega, Nm1, 1, G[1], G[0]);
			return begin + 2; // removing first iteration from it
		}
		return begin;
	}
	public final void call_barrier(final int id, final int total_threads) {
		// Signal this thread has done iteration
		sync[id][0]++;

		// Wait for neighbors
		if (id > 0) {
			while (sync[id-1][0] < sync[id][0]) ;
		}
		if (id < total_threads -1) {
			while (sync[id+1][0] < sync[id][0]) ;
		}
	}

	private void middle_rows(final double G[][], final double omega_4, final double one_m_omega, 
							final int Nm1, int i, final double[] Gi, final double[] Gim1) 
	{
		final double [] Gip1 = G[i+1];
		final double [] Gim2 = G[i-2];
		int j;
		final int Nm3 = Nm1-2;
		
		// removing the if((j+1) != Nm1)  condition from the inside loop body
		for (j = 1; j < Nm3; j += 2)
		{
			Gi[j] 	  = omega_4 * (Gim1[j] + Gip1[j] + Gi[j-1] + Gi[j+1]) + one_m_omega * Gi[j];
			Gim1[j+1] = omega_4 * (Gim2[j+1] + Gi[j+1] + Gim1[j] + Gim1[j+2]) + one_m_omega * Gim1[j+1];
		}
		
		Gi[j] = omega_4 * (Gim1[j] + Gip1[j] + Gi[j-1] + Gi[j+1]) + one_m_omega * Gi[j];
		
		if((j+1) != Nm1) 
		{
			Gim1[j+1] = omega_4 * (Gim2[j+1] + Gi[j+1] + Gim1[j]+ Gim1[j+2]) + one_m_omega * Gim1[j+1];
		}
	}

	private void last_row(	final double G[][], final double omega_4, final double one_m_omega, 
							final int Nm1, int i, final double[] Gi, final double[] Gim1) {
		final double [] Gim2 = G[i-2];

		// removing if((j+1) != Nm1) condition using Nm-1 instead of Nm1
		final int Nm2 = Nm1 - 1;
		
		for (int j = 1; j < Nm2; j += 2)
		{
			Gim1[j+1] = omega_4 * (Gim2[j+1] + Gi[j+1] + Gim1[j] + Gim1[j+2]) + one_m_omega * Gim1[j+1];
		}
	}

	private static void first_row (	final double G[][], 	  final double omega_over4, final double one_m_omega, 
			final int Nm1, int i, final double[] Gi, final double[] Gim1) 
	{
		final double [] Gip1 = G[i+1];

		for (int j = 1; j < Nm1; j += 2)
		{
			Gi[j] = omega_over4 * (Gim1[j] + Gip1[j] + Gi[j-1] + Gi[j+1]) + one_m_omega * Gi[j];
		}
	}
}
