//#include <iostream>
//#include <fstream>
//#include <string>
//#include <mpi.h>
//#include <omp.h>
//#include <stdio.h>
//#include <stdlib.h>
//using namespace std;
//int NUM_OF_TEST_CASES = 4;
//int main(int argc, char* argv[])
//{
//
//	string arr[4] = { "C://Users/asd/source/repos/Project1/test case1.txt",
//		"C://Users/asd/source/repos/Project1/test case2.txt",
//		"C://Users/asd/source/repos/Project1/test case3 .txt",
//		"C://Users/asd/source/repos/Project1/test case4 .txt"
//	};
//	string substring = "GCCAGATATTCCCCCCGTT";
//	int count = 0;
//	int rank, size;
//	int choice = 100;
//	string line;
//	ifstream Myfile;
//
//
//	MPI_Init(&argc, &argv);
//	MPI_Comm_size(MPI_COMM_WORLD, &size);
//	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
//
//	if (rank == 0)
//	{
//		cout << "Parallel Task!\n";
//		cout << "[1] To run first test case\n";
//		cout << "[2] To run second test case\n";
//		cout << "[3] To run third test case\n";
//		cout << "[4] To run fourth test case\n";
//
//		cin >> choice;
//		if (choice <= 4)
//		{
//			Myfile.open(arr[choice - 1]);
//			if (Myfile.is_open())
//			{
//				getline(Myfile, line);
//
//				Myfile.close();
//			}
//			else cout << "ERROR" << endl;
//			MPI_Send(line.c_str(), line.size(), MPI_CHAR, 1, 0, MPI_COMM_WORLD);
//
//		}
//
//	}
//	else if (rank == 1)
//	{
//		MPI_Status x_Status;
//		int count;
//		int j, k, id;
//		int subcount = 0;
//		int countA = 0;
//		int countC = 0;
//		int countG = 0;
//		int countT = 0;
// 
//		MPI_Probe(0, 0, MPI_COMM_WORLD, &x_Status);
//		MPI_Get_count(&x_Status, MPI_CHAR, &count);
//		char* buffer = new char[count];
//		MPI_Recv(buffer, count, MPI_CHAR, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
//		omp_set_num_threads(5);
//		int Chunksize = count / 5;
//#pragma omp parallel shared (countA,countC,countG,countT,subcount,substring) private(id,j,k)
//		{
//			id = omp_get_thread_num();
//			int start = id * Chunksize;
//			int end = start + Chunksize;
//          long long microseconds_seq;
//          long long microseconds_char;
//			int local_subcount = 0;
//			int localsubcount = 0;
//			int local_counterA = 0;
//			int local_counterC = 0;
//			int local_counterG = 0;
//			int local_counterT = 0;
//
//			for (j = start; j < end; j++)
//			{
//              auto start_char = std::chrono::high_resolution_clock::now();
// 	            
//				if (buffer[j] == 'A')
//				{
//					local_counterA++;
//				}
//				else if (buffer[j] == 'C')
//				{
//					local_counterC++;
//				}
//				else if (buffer[j] == 'G')
//				{
//					local_counterG++;
//				}
//				else if (buffer[j] == 'T')
//				{
//					local_counterT++;
//				}
				
//              auto start_seq = std::chrono::high_resolution_clock::now();
//				if (buffer[j] == substring[0])
//				{
//					bool found = true;
//					for (k = 1; k < substring.size(); k++)
//					{
//						if (k + j >= count || buffer[k + j] != substring[k])
//						{
//							found = false;
//							break;
//						}
//					}
//
//					if (found)
//					{
//						local_subcount++;
//					}
//				}
//				auto elapsed_seq = std::chrono::high_resolution_clock::now() - start_seq;
//			microseconds_seq = std::chrono::duration_cast<std::chrono::microseconds>(
// 				elapsed_seq).count(); 
//			}
//#pragma omp critical
//			{
//				countA += local_counterA;
//				countC += local_counterC;
//				countG += local_counterG;
//				countT += local_counterT;
//				subcount += local_subcount;
//
//			}
//
//		}
//		cout << "A appeared " << countA << " Times" << endl;
//		cout << "C appeared " << countC << " Times" << endl;
//		cout << "G appeared " << countG << " Times" << endl;
//		cout << "T appeared " << countT << " Times" << endl;
//		cout << substring << " appeared " << subcount << " times" << endl;
//      cout << "Test Case"<<arr[0]<< Finished\nTime for Characters: " << microseconds_char <<endl<<"Time for Sequence: "<<microseconds_seq << endl;
//	}
//	MPI_Finalize();
//	return 0;
//
//}

//#include <iostream>
//#include <omp.h>
//#include <stdio.h>
//#include <stdlib.h>
//#include <vector>
//#include <fstream>
//#include <string>
//using namespace std;
//int main()
//{
//
//	string substring = "GCCAGAT";
//	int j, k, id;
//	int count = 0;
//	int countA = 0;
//	int countC = 0;
//	int countG = 0;
//	int countT = 0;
//	int subcount = 0;
//	ifstream Myfile;
//	char* buffer = new char[count];
//	string arr[4] = { 
//		"C://Users/asd/source/repos/Project1/test case1.txt",
//		"C://Users/asd/source/repos/Project1/test case2.txt",
//		"C://Users/asd/source/repos/Project1/test case3 .txt",
//		"C://Users/asd/source/repos/Project1/test case4 .txt"
//	};
//	string line;
//	Myfile.open(arr[3]);
//	if (Myfile.is_open())
//	{
//		getline(Myfile, line);
//
//		Myfile.close();
//	}
//	else cout << "Unable to open this file." << endl;
//	omp_set_num_threads(5);
//	int ChunkSize = line.size() / 5;
//	#pragma omp parallel shared (countA,countC,countG,countT,subcount,substring) private(id,j,k)
//	{
//		id = omp_get_thread_num();
//		int start = id * ChunkSize;
//		int end = start + ChunkSize;
//		int localsubcount = 0;
//		int local_counterA = 0;
//		int local_counterC = 0;
//		int local_counterG = 0;
//		int local_counterT = 0;
//
//		for (j = start; j <= end - 1; j++)
//		{
//
//			if (line[j] == 'A')
//			{
//				local_counterA++;
//			}
//			else if (line[j] == 'C')
//			{
//				local_counterC++;
//			}
//			else if (line[j] == 'G')
//			{
//				local_counterG++;
//			}
//			else if (line[j] == 'T')
//			{
//				local_counterT++;
//			}
//			if (line[j] == substring[0])
//			{
//				bool found = true;
//				for (k = 1; k < substring.size(); k++)
//				{
//					
//					if (k + j >= line.size() || line[k + j] != substring[k])
//					{
//						found = false;
//						break;
//					}
//				}
//				if (found)
//				{
//					localsubcount++;
//				}
//			}
//		}
//		#pragma omp critical
//		{
//			countA += local_counterA;
//			countC += local_counterC;
//			countG += local_counterG;
//			countT += local_counterT;
//			subcount += localsubcount;
//
//		}
//
//
//
//	}
//	cout << "A appeared " << countA << " Times" << endl;
//	cout << "C appeared " << countC << " Times" << endl;
//	cout << "G appeared " << countG << " Times" << endl;
//	cout << "T appeared " << countT << " Times" << endl;
//
//	cout << substring << " appeared " << subcount << " times" << endl;
//
//	return 0;
//}

//#include <iostream>
//#include <fstream>
//#include <string>
//#include <mpi.h>
//#include<chrono>
//using namespace std;
//int NUM_OF_TEST_CASES = 4;
//int main(int argc, char* argv[])
//{
//
//	string arr[4] = { "C://Users/asd/source/repos/Project1/test case1.txt",
//		"C://Users/asd/source/repos/Project1/test case2.txt",
//		"C://Users/asd/source/repos/Project1/test case3 .txt",
//		"C://Users/asd/source/repos/Project1/test case4 .txt"
//	};
//	string substring = "GCCAGATATTCCCCCCGTT";
//	int count = 0;
//	int rank, size;
//	int choice = 100;
//	string line;
//	ifstream Myfile;
//	
//
//
//	MPI_Init(&argc, &argv);
//	
//	MPI_Comm_size(MPI_COMM_WORLD, &size);
//	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
//
//	if (rank == 0)
//	{
//		cout << "Parallel Task!\n";
//		cout << "[1] To run first test case\n";
//		cout << "[2] To run second test case\n";
//		cout << "[3] To run third test case\n";
//		cout << "[4] To run fourth test case\n";
//		cout << "[5] To run all test cases\n";
//		cin >> choice;
//		if (choice == 5)
//		{
//			for (int i = 0; i < NUM_OF_TEST_CASES; i++)
//			{
//				Myfile.open(arr[i]);
//				if (Myfile.is_open())
//				{
//					getline(Myfile, line);
//
//					Myfile.close();
//				}
//				else cout << "ERROR" << endl;
//				MPI_Send(line.c_str(), line.size() + 1, MPI_CHAR, i + 1, 0, MPI_COMM_WORLD);
//			}
//		}
//		else {
//			Myfile.open(arr[choice - 1]);
//			if (Myfile.is_open())
//			{
//				getline(Myfile, line);
//
//				Myfile.close();
//			}
//			else cout << "ERROR" << endl;
//			MPI_Send(line.c_str(), line.size() + 1, MPI_CHAR, 1, 0, MPI_COMM_WORLD);
//			for (int i = 1; i < NUM_OF_TEST_CASES; i++)
//			{
//				string x = "x";
//				MPI_Send(x.c_str(), 1, MPI_CHAR, i + 1, 0, MPI_COMM_WORLD);
//
//			}
//		}
//	}
//	else if (rank == 1)
//	{
//		auto start = std::chrono::high_resolution_clock::now();
//		MPI_Status x_Status;
//		int count;
//		int countA = 0;
//		int countC = 0;
//		int countG = 0;
//		int countT = 0;
//		MPI_Probe(0, 0, MPI_COMM_WORLD, &x_Status);
//		MPI_Get_count(&x_Status, MPI_CHAR, &count);
//		char* buffer = new char[count];
//		MPI_Recv(buffer, count, MPI_CHAR, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
//
//		cout << "Recived!" << endl;
//
//		for (int j = 0; j < count; j++)
//		{
//
//			if (buffer[j] == 'A')
//			{
//				countA++;
//			}
//			else if (buffer[j] == 'C')
//			{
//				countC++;
//			}
//			else if (buffer[j] == 'G')
//			{
//				countG++;
//			}
//			else if (buffer[j] == 'T')
//			{
//				countT++;
//			}
//		}
//
//
//		cout << "A appeared " << countA << " Times" << endl;
//		cout << "C appeared " << countC << " Times" << endl;
//		cout << "G appeared " << countG << " Times" << endl;
//		cout << "T appeared " << countT << " Times" << endl;
//
//
//		int subcount = 0;
//		for (int j = 0; j < count; j++)
//		{
//			if (buffer[j] == substring[0])
//			{
//				bool found = true;
//				for (int k = 1; k < substring.size(); k++)
//				{
//					if (k + j >= count || buffer[k + j] != substring[k])
//					{
//						found = false;
//						break;
//					}
//				}
//
//				if (found)
//				{
//					subcount++;
//				}
//			}
//		}
//
//		cout << substring << " appeared " << subcount << " times" << endl;
//		double stop_SUB = MPI_Wtime();
//		double TotalTime_SUB = (stop_SUB - start_sub) / 1000;
//		cout << "TIME TO COMPUTE SEQUENCE: " << TotalTime_SUB << endl;
//	}
//	else if (rank == 2)
//	{
//		MPI_Status x_Status;
//		int count;
//		int countA = 0;
//		int countC = 0;
//		int countG = 0;
//		int countT = 0;
//		MPI_Probe(0, 0, MPI_COMM_WORLD, &x_Status);
//		MPI_Get_count(&x_Status, MPI_CHAR, &count);
//		char* buffer = new char[count];
//		MPI_Recv(buffer, count, MPI_CHAR, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
//		if (count != 1)
//		{
//			cout << "Recived!" << endl;
//
//			for (int j = 0; j < count; j++)
//			{
//
//				if (buffer[j] == 'A')
//				{
//					countA++;
//				}
//				else if (buffer[j] == 'C')
//				{
//					countC++;
//				}
//				else if (buffer[j] == 'G')
//				{
//					countG++;
//				}
//				else if (buffer[j] == 'T')
//				{
//					countT++;
//				}
//			}
//
//
//			cout << "A appeared " << countA << " Times" << endl;
//			cout << "C appeared " << countC << " Times" << endl;
//			cout << "G appeared " << countG << " Times" << endl;
//			cout << "T appeared " << countT << " Times" << endl;
//
//
//
//			int subcount = 0;
//			for (int j = 0; j < count; j++)
//			{
//				if (buffer[j] == substring[0])
//				{
//					bool found = true;
//					for (int k = 1; k < substring.size(); k++)
//					{
//						if (k + j >= count || buffer[k + j] != substring[k])
//						{
//							found = false;
//							break;
//						}
//					}
//
//					if (found)
//					{
//						subcount++;
//					}
//				}
//			}
//
//			cout << substring << " appeared " << subcount << " times" << endl;
//
//		}
//	}
//	else if (rank == 3)
//	{
//		MPI_Status x_Status;
//		int count;
//		int countA = 0;
//		int countC = 0;
//		int countG = 0;
//		int countT = 0;
//		MPI_Probe(0, 0, MPI_COMM_WORLD, &x_Status);
//		MPI_Get_count(&x_Status, MPI_CHAR, &count);
//		char* buffer = new char[count];
//		MPI_Recv(buffer, count, MPI_CHAR, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
//		if (count != 1)
//		{
//			cout << "Recived!" << endl;
//
//			for (int j = 0; j < count; j++)
//			{
//
//				if (buffer[j] == 'A')
//				{
//					countA++;
//				}
//				else if (buffer[j] == 'C')
//				{
//					countC++;
//				}
//				else if (buffer[j] == 'G')
//				{
//					countG++;
//				}
//				else if (buffer[j] == 'T')
//				{
//					countT++;
//				}
//			}
//
//
//			cout << "A appeared " << countA << " Times" << endl;
//			cout << "C appeared " << countC << " Times" << endl;
//			cout << "G appeared " << countG << " Times" << endl;
//			cout << "T appeared " << countT << " Times" << endl;
//
//
//
//			int subcount = 0;
//			for (int j = 0; j < count; j++)
//			{
//				if (buffer[j] == substring[0])
//				{
//					bool found = true;
//					for (int k = 1; k < substring.size(); k++)
//					{
//						if (k + j >= count || buffer[k + j] != substring[k])
//						{
//							found = false;
//							break;
//						}
//					}
//
//					if (found)
//					{
//						subcount++;
//					}
//				}
//			}
//
//			cout << substring << " appeared " << subcount << " times" << endl;
//		}
//	}
//	else if (rank == 4)
//	{
//		MPI_Status x_Status;
//		int count;
//		int countA = 0;
//		int countC = 0;
//		int countG = 0;
//		int countT = 0;
//		MPI_Probe(0, 0, MPI_COMM_WORLD, &x_Status);
//		MPI_Get_count(&x_Status, MPI_CHAR, &count);
//		char* buffer = new char[count];
//		MPI_Recv(buffer, count, MPI_CHAR, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
//		if (count != 1)
//		{
//			cout << "Recived!" << endl;
//
//			for (int j = 0; j < count; j++)
//			{
//
//				if (buffer[j] == 'A')
//				{
//					countA++;
//				}
//				else if (buffer[j] == 'C')
//				{
//					countC++;
//				}
//				else if (buffer[j] == 'G')
//				{
//					countG++;
//				}
//				else if (buffer[j] == 'T')
//				{
//					countT++;
//				}
//			}
//
//
//			cout << "A appeared " << countA << " Times" << endl;
//			cout << "C appeared " << countC << " Times" << endl;
//			cout << "G appeared " << countG << " Times" << endl;
//			cout << "T appeared " << countT << " Times" << endl;
//			
//
//			double start_sub = MPI_Wtime();
//			int subcount = 0;
//			for (int j = 0; j < count; j++)
//			{
//				if (buffer[j] == substring[0])
//				{
//					bool found = true;
//					for (int k = 1; k < substring.size(); k++)
//					{
//						if (k + j >= count || buffer[k + j] != substring[k])
//						{
//							found = false;
//							break;
//						}
//					}
//
//					if (found)
//					{
//						subcount++;
//					}
//				}
//			}
//			
//			cout << substring << " appeared " << subcount << " times" << endl;
//		
//		}
//	}
//	
//	MPI_Finalize();
//
//	return 0;
//
//}

//#include <iostream>
//#include <fstream>
//#include <string>
//#include<vector>
//#include<chrono>
//
//using namespace std;
//vector<string> readFile(string filename)
//{
//	ifstream Myfile;
//	vector<string> lines;
//	Myfile.open(filename);
//
//	if (Myfile.is_open())
//	{
//		string line;
//
//		while (getline(Myfile, line))
//			lines.push_back(line);
//
//		Myfile.close();
//	}
//	else
//		cout << "Unable to open this file." << endl;
//	return lines;
//}
//long long microseconds_char;
//void countACGT(string filename)
//{
//	//double start, Stop, TotalTime = 0.0;
//	//start = clock();
//
//	vector<string> file = readFile(filename);
//	int countA = 0;
//	int countC = 0;
//	int countG = 0;
//	int countT = 0;
//
//	auto start_char = std::chrono::high_resolution_clock::now();
//
//
//	for (string line : file)
//	{
//		for (char charachter : line)
//		{
//
//			if (charachter == 'A')
//			{
//
//				countA++;
//
//			}
//			else if (charachter == 'C')
//			{
//				countC++;
//			}
//			else if (charachter == 'G')
//			{
//				countG++;
//			}
//			else if (charachter == 'T')
//			{
//				countT++;
//			}
//		}
//	}
//
//
//	cout << "A appeared " << countA << " Times" << endl;
//	cout << "C appeared " << countC << " Times" << endl;
//	cout << "G appeared " << countG << " Times" << endl;
//	cout << "T appeared " << countT << " Times" << endl;
//
//	auto elapsed_char = std::chrono::high_resolution_clock::now() - start_char;
//	microseconds_char = std::chrono::duration_cast<std::chrono::microseconds>(
//		elapsed_char).count();
//}
//long long microseconds_seq;
//void countOccurence(string filename, string substring)
//{
//	//double start_SUB, Stop_SUB = 0.0;
//	//double TotalTime_SUB = 0.0;
//
//
//
//	auto start_seq = std::chrono::high_resolution_clock::now();
//	vector<string> file = readFile(filename);
//	string::size_type pos = 0;
//	int count = 0;
//	for (string line : file) {
//		while ((pos = line.find(substring, pos)) != string::npos)
//		{
//			count++;
//			pos += substring.length();
//		}
//	}
//
//
//	cout << substring << " appeared " << count << " times" << endl;
//	auto elapsed_seq = std::chrono::high_resolution_clock::now() - start_seq;
//	microseconds_seq = std::chrono::duration_cast<std::chrono::microseconds>(
//		elapsed_seq).count();
//
//}
//int main()
//{
//
//
//
//	/*cout << "[Clock] Started" << endl << endl;*/
//
//	cout << "\t \t \t \t \t" << "##########################" << endl;
//	cout << "\t \t \t \t \t" << "##" << "\t" << "test case 1" << "\t" << "##" << endl;
//	cout << "\t \t \t \t \t" << "##########################" << endl;
//	countACGT("test case1.txt");
//	cout << endl;
//	countOccurence("test case1.txt", "GCCAGATATTCCCCCCGTT");
//	cout << endl;
//
//
//	//check mess, ok sec
//
//	cout << "Test Case 1 Finished\nTime for Characters: " << microseconds_char <<endl<<"Time for Sequence: "<<microseconds_seq << endl;
//
//
//	cout << "\t \t \t \t \t" << "##########################" << endl;
//	cout << "\t \t \t \t \t" << "##" << "\t" << "test case 2" << "\t" << "##" << endl;
//	cout << "\t \t \t \t \t" << "##########################" << endl;
//	countACGT("test case2.txt");
//	cout << endl;
//	countOccurence("test case2.txt", "GCCAGATATTCCCCCCGTT");
//	cout << endl;
//
//	cout << "Test Case 2 Finished\nTime for Characters: " << microseconds_char << endl << "Time for Sequence: " << microseconds_seq << endl;
//	
//
//
//
//
//	cout << "\t \t \t \t \t" << "##########################" << endl;
//	cout << "\t \t \t \t \t" << "##" << "\t" << "test case 3" << "\t" << "##" << endl;
//	cout << "\t \t \t \t \t" << "##########################" << endl;
//	countACGT("test case3 .txt");
//	cout << endl;
//	countOccurence("test case3 .txt", "GCCAGATATTCCCCCCGTT");
//	cout << endl;
//
//	cout << "Test Case 3 Finished\nTime for Characters: " << microseconds_char << endl << "Time for Sequence: " << microseconds_seq << endl;
//	
//
//	
//
//	cout << "\t \t \t \t \t" << "##########################" << endl;
//	cout << "\t \t \t \t \t" << "##" << "\t" << "test case 4" << "\t" << "##" << endl;
//	cout << "\t \t \t \t \t" << "##########################" << endl;
//	countACGT("test case4 .txt");
//	cout << endl;
//	countOccurence("test case4 .txt", "GCCAGATATTCCCCCCGTT");
//	cout << endl;
//
//	
//
//	cout << "Test Case 4 Finished\nTime for Characters: " << microseconds_char << endl << "Time for Sequence: " << microseconds_seq << endl;
//
//
//	return 0;
//}
