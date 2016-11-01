#include<bits/stdc++.h>
#include <iostream>
#include <cmath>
#include <set>
#include "UCTGraph.h"
double p, q, r;
int timeStamp;
vector< set < pair<int, pair<int, int>  > > > lift1Stopped;   // 1: floor number 2.1:time at which button up was presssed  2.2 time at which button down was prssed  
vector< set < pair<int, pair<int, int>  > > > lift2Stopped;
vector<pair<State, double>> nextStates ;
map<ElevatorAction,string> decodes;
double P[
		T + 2][
		N + 1][
		N + 1][
		T + 2][T + 2];  // prob. of n people going from f1 to f2 after t1 - up waiting time, t2 - down waiting time
double P_arrived[N + 1][T + 2][T + 2];
double n_exp[
		N + 1][
		N + 1][
		T + 2][
		T + 2];   // expected no. of people going from f1 to f2 after t1 - up waiting time, t2 - down waiting time
double n_exp_up[N + 1][T + 2][T + 2];     // expected no. of people going up from f
double n_exp_down[N + 1][T + 2][T + 2];   // expected no. of people going down from f

int d[N + 1][N + 1][2][2];  // distance between elevator and person given direction of elevator and person


void precompute_P_arrival(double p, double q, double r) {
	memset(P_arrived, 0, sizeof(P_arrived));

	for (int floor = 1; floor <= N; ++floor) {
		P_arrived[floor][0][0] = 1;     // time passed = 0
		for (int t = 1; t <= T; ++t) {  // time t
			for (int m = 0; m <= t; ++m) {
				double prob = 0;
				double arrival_choice_factor;
				if (floor == 1) {
					arrival_choice_factor = q;
				}
				else {
					arrival_choice_factor = (1.0 - q) / (double) N;
				}
				for (int k = m; k <= t; ++k) {
					prob += ((double) C[t][k] * pow(p, k) * pow((1.0 - p), (t - k)) * (double) C[k][m] *
					         pow(arrival_choice_factor, m) * pow((1.0 - arrival_choice_factor), (k - m)));
				}
				P_arrived[floor][m][t] = prob;
			}
		}
	}
}

void precompute_P(double p, double q, double r) {
	memset(P, 0, sizeof(P));
	for (int t1 = 0; t1 <= T; ++t1) {
		for (int t2 = 0; t2 <= T; ++t2) {
			for (int floor_start = 1; floor_start <= N; ++floor_start) {
				// for floor_end = floor_start, probabilities are 0
				// For t = 0, all probabilities are 0
				for (int floor_end = 1; floor_end <= N; ++floor_end) {
					if (floor_end == floor_start) {
						continue;
					}

					double prob;
					double alight_choice_factor;
					if (floor_start == 1) {
						if (t1 != 0) {
							alight_choice_factor = 1.0 / (double) (N - 1);
						}
						else {
							alight_choice_factor = 0.0;
						}
					}
					else {
						if (t1 == 0 and t2 == 0) {
							alight_choice_factor = 0.0;
						}
						else if (t1 == 0 and t2 != 0) {
							if (floor_end > floor_start) {
								alight_choice_factor = 0.0;
							}
							else if (floor_end < floor_start) {
								if (floor_end == 1) {
									alight_choice_factor =
											r / (r + (double) (floor_start - 1) * ((1.0 - r) / (N - 2.0)));
								}
								else {
									alight_choice_factor = ((1.0 - r) / (N - 2.0)) /
									                       (r + (double) (floor_start - 1) * ((1.0 - r) / (N - 2.0)));
								}
							}
						}
						else if (t1 != 0 and t2 == 0) {
							if (floor_end > floor_start) {
								alight_choice_factor = 1.0 / (double) (N - floor_start);
							}
							else if (floor_end < floor_start) {
								alight_choice_factor = 0.0;
							}
						}
						else if (t1 != 0 and t2 != 0) {
							if (floor_end == 1) {
								alight_choice_factor = r;
							}
							else {
								alight_choice_factor = (1.0 - r) / (N - 2.0);
							}
						}
					}

					if (floor_end > floor_start) {   // Up
						for (int n = 0; n <= t1; ++n) {
							prob = 0;
							for (int m = max(1, n); m <= (t1 - 1) + 1; ++m) {
								prob += ((double) C[m][n] * pow(alight_choice_factor, n) *
								         pow(1.0 - alight_choice_factor, (m - n)) *
								         P_arrived[floor_start][m - 1][max(t1 - 1, 0)]);
							}
							P[n][floor_start][floor_end][t1][t2] = prob;
						}
					}
					else if (floor_end < floor_start) {  // Down
						for (int n = 0; n <= t2; ++n) {
							prob = 0;
							for (int m = max(1, n); m <= (t2 - 1) + 1; ++m) {
								prob += ((double) C[m][n] * pow(alight_choice_factor, n) *
								         pow(1.0 - alight_choice_factor, (m - n)) *
								         P_arrived[floor_start][m - 1][max(t2 - 1, 0)]);
							}
							P[n][floor_start][floor_end][t1][t2] = prob;
						}
					}
				}
			}
		}
	}
}

void precompute_n_exp(double p, double q, double r) {
	memset(n_exp, 0, sizeof(n_exp));

	for (int t1 = 0; t1 <= T; ++t1) {
		for (int t2 = 0; t2 <= T; ++t2) {
			for (int floor_start = 1; floor_start <= N; ++floor_start) {
				// for floor_end = floor_start, n_exp = 0
				for (int floor_end = 1; floor_end <= N; ++floor_end) {
					if (floor_start == floor_end) {
						continue;
					}
					if (floor_end > floor_start) {
						double exp = 0;
						for (int n = 1; n <= (t1 - 1) + 1; ++n) {
							exp += ((double) (n) * P[n][floor_start][floor_end][t1][t2]);
						}
						n_exp[floor_start][floor_end][t1][t2] = exp;
					}
					else if (floor_end < floor_start) {
						double exp = 0;
						for (int n = 1; n <= (t2 - 1) + 1; ++n) {
							exp += ((double) (n) * P[n][floor_start][floor_end][t1][t2]);
						}
						n_exp[floor_start][floor_end][t1][t2] = exp;
					}
				}
			}
		}
	}
}

void precompute_n_exp_up_down() {
	memset(n_exp_up, 0, sizeof(n_exp_up));
	memset(n_exp_down, 0, sizeof(n_exp_down));

	for (int floor_start = 1; floor_start <= N; ++floor_start) {
		for (int t1 = 0; t1 <= T; ++t1) {
			for (int t2 = 0; t2 <= T; ++t2) {
				// Down computation
				double exp_down = 0;
				for (int floor_end = 1; floor_end < floor_start; ++floor_end) {
					exp_down += n_exp[floor_start][floor_end][t1][t2];
				}
				n_exp_down[floor_start][t1][t2] = exp_down;

				// Up computation
				double exp_up = 0;
				for (int floor_end = floor_start + 1; floor_end <= N; ++floor_end) {
					exp_up += n_exp[floor_start][floor_end][t1][t2];
				}
				n_exp_up[floor_start][t1][t2] = exp_up;
			}
		}
	}
}

void precompute_distance() {
	memset(d, 0, sizeof(d));
	for (int i = 1; i <= N; ++i) {  // elevator position
		// elevator up, person up
		for (int j = 1; j < i; ++j) {
			d[i][j][1][1] = N - i + N - 1 + j - 1;
		}
		for (int j = i; j <= N; ++j) {
			d[i][j][1][1] = (j - i);
		}

		// elevator up, person down
		for (int j = 1; j <= N; ++j) {
			d[i][j][1][0] = N - i + N - j;
		}

		// elevator down, person up
		for (int j = 1; j <= N; ++j) {
			d[i][j][0][1] = i - 1 + j - 1;
		}

		// elevator down, person down
		for (int j = 1; j <= i; ++j) {
			d[i][j][0][0] = i - j;
		}
		for (int j = i + 1; j <= N; ++j) {
			d[i][j][0][0] = i - 1 + N - 1 + N - j;
		}
	}
}

void precompute(double p, double q, double r) {
	precompute_P_arrival(p, q, r);
	precompute_P(p, q, r);
	precompute_n_exp(p, q, r);
	precompute_n_exp_up_down();
	precompute_distance();
	for(int i=1;i<N+1;i++)
	{
		set<pair<int , pair< int , int>  >> temp ;
		lift1Stopped.push_back(temp) ;
		lift2Stopped.push_back(temp) ;
 	}
 	decodes[AU]="AU";
  	decodes[AD]="AD";
  	decodes[AOU]="AOU";
  	decodes[AOD]="AOD";
}

double getRemCost(Elevator elevator) {
	double rem_cost = 0.0;
	if (elevator.is_up and elevator.elevatorState != EMPTY) {
		for (int i = elevator.position + 1; i <= N; ++i) {
			rem_cost += ((double) (i - elevator.position) * elevator.alight[i]);
		}
	}
	else if (!elevator.is_up and elevator.elevatorState != EMPTY) {
		for (int i = elevator.position - 1; i >= 1; --i) {
			rem_cost += (double) ((elevator.position - i) * elevator.alight[i]);
		}
	}
	return rem_cost;
}

double evaluateTerminalState(State state) {
	// cost of remaining people in elevator
	double rem_cost = 0;
	rem_cost += getRemCost(state.elevator1);
	rem_cost += getRemCost(state.elevator2);

	// cost of waiting people
	double wait_cost = 0;
	for (int i = 1; i <= N; ++i) {

	}

	return (rem_cost + wait_cost);
}

//// returns cost of performing given action sequences on elevators in current state
//double getSequenceCost(State state, vector<ElevatorAction> action_seq1, vector<ElevatorAction> action_seq2) {
//	assert(action_seq1.size() == action_seq2.size());   // action sequences of same size
//	int num_actions = action_seq1.size();
//	for (int i = 0; i < num_actions; ++i) {
//
//	}
//}
 
 void findNextStates2(int floor1, int floor2, State s, double prob)
 {

 	if(floor2>=N+1)
 		{	nextStates.push_back(pair<State , double>(s, prob));
 			return; 
 		}
 	if(floor1>=N+1)
 		{
 			if(s.elevator2.btnPressed[floor2])
 			{
 				for(int i=1;i<N+1;i++)
	 			{
	 				State temp= s; 
	 				temp.elevator2.alight[floor2]=i ;
	 				double p=0.0 ;
	 				set<pair<int, pair<int ,int>>>::iterator it ;
	 				for (it=lift2Stopped[floor2].begin(); it!=lift2Stopped[floor2].end()&&lift2Stopped[floor2].size(); ++it)
    				{//cout << ' ' << it->first<<it->second;
    					p= p+ P[i][it->first][floor2][it->second.first][it->second.second] ;
    				}	
	 				findNextStates2(floor1, floor2+1, temp, prob*p);
	 			}	
 			}
 			else
 			{
 				findNextStates2(floor1,floor2+1,s,prob);		
 			}
 		}	
 	else if(s.elevator1.btnPressed[floor1])
 		{
 			for(int i=1;i<N+1;i++)
 			{	 
 				State temp= s; 
 				temp.elevator1.alight[floor1]=i ;
 				
 				double p=0.0 ;
	 			set<pair<int, pair<int ,int>>>::iterator it ;
	 			for (it=lift1Stopped[floor1].begin(); it!=lift1Stopped[floor1].end()&&lift1Stopped[floor1].size(); ++it)
    			{//cout << ' ' << it->first<<it->second;
    				cout<<"yo"<<endl;
    				p= p+ P[i][it->first][floor2][it->second.first][it->second.second] ;
    			}
 				findNextStates2(floor1+1, floor2, temp, prob*p);
 			}
 		}
 		else
 		{
 			findNextStates2(floor1+1,floor2,s,prob);
 		}	 
 } 

vector< pair<Action,double>> dummyCost()
{
	vector< pair<Action,double>> temp ;
	for(int i=0;i<4;i++)
	{
		pair<Action,double> t1=pair<Action,double>( pair<ElevatorAction,ElevatorAction>(AU,AU),100); 	
		temp.push_back(t1) ;
	}
	return temp ;
}
void changes(State &s, ElevatorAction s1, ElevatorAction s2)
{
	for(int i=1;i<N+1;i++)
	{
		if(s.time_up[i]>0)
		{
			s.time_up[i]++; 
		}
		if(s.time_down[i]>0)
		{
			s.time_down[i]++; 
		}
	}
	switch(s1)
	{
		case AU :
		break;
		case AOU : lift1Stopped[s.elevator1.position].clear();  
					for(int i=s.elevator1.position+1;i<N+1;i++)  
						{
							int a=s.time_up[s.elevator1.position] ;
							int b=s.time_down[s.elevator1.position]; 
							pair<int, int > temp= pair<int, int >(a,b) ;
							lift1Stopped[i].insert(pair<int, pair<int , int> >(s.elevator1.position,temp)); 
			
						}
						s.time_up[s.elevator1.position]=0;
						s.elevator1.btnPressed[s.elevator1.position]=false; 
		break ;
		case AD :
		break; 
		case AOD :  lift1Stopped[s.elevator1.position].clear();  
					for(int i=s.elevator1.position-1;i>0;i--)  
						{
							int a=s.time_up[s.elevator1.position] ;
							int b=s.time_down[s.elevator1.position]; 
							pair<int, int > temp= pair<int, int >(a,b) ;
							lift1Stopped[i].insert(pair<int, pair<int , int> >(s.elevator1.position,temp)); 
			
						}
						s.time_down[s.elevator1.position]=0;
						s.elevator1.btnPressed[s.elevator1.position]=false;
		break;
		default:  break ;

	}
	switch(s2)
	{
		case AU :
		break;
		case AOU :lift2Stopped[s.elevator2.position].clear();  
					for(int i=s.elevator2.position+1;i<N+1;i++)  
						{
							int a=s.time_up[s.elevator2.position] ;
							int b=s.time_down[s.elevator2.position]; 
							pair<int, int > temp= pair<int, int >(a,b) ;
							lift2Stopped[i].insert(pair<int, pair<int , int> >(s.elevator2.position,temp)); 
			
						}
						s.time_up[s.elevator2.position]=0;
						s.elevator2.btnPressed[s.elevator2.position]=false;
		break ;
		case AD :
		break; 
		case AOD : lift2Stopped[s.elevator2.position].clear();  
					for(int i=s.elevator2.position-1;i>0;i--)  
						{
							int a=s.time_up[s.elevator2.position] ;
							int b=s.time_down[s.elevator2.position]; 
							pair<int, int > temp= pair<int, int >(a,b) ;
							lift2Stopped[i].insert(pair<int, pair<int , int> >(s.elevator2.position,temp)); 
			
						}
						s.time_down[s.elevator2.position]=0;
						s.elevator2.btnPressed[s.elevator2.position]=false;
		break;
		default: break;

	}
} 
Action findAction()
{
	double expectiCost[25];   // 0 : AU   1:AD 
	memset(expectiCost, 0, sizeof(expectiCost));
	vector<pair<Action,double> > v;
	for(int i=0;i<nextStates.size();i++)
	{	State *s123=  & nextStates[i].first;
		//printS(s123) ;
		v=UCTGraph::getBaseCosts(s123,1);   //dummyCost();//call fxn nextStates[i].second* //
		for(int j=0;j<v.size();j++)
		{
			expectiCost[j] = expectiCost[j]+ nextStates[i].second*v[j].second ;
		} 
	} 
	double temp = LONG_MIN ;
	int act ;
	for(int i=0;i<v.size();i++)
	{
		if(expectiCost[i]>temp)
		{
			act= i ; 
		}
	} 
	return v[act].first ;
}
int main() {
	p = 0.8;
	q = 0.5;
	r = 0.5;
	precompute(p, q, r);
   cerr<<"started";
//	 find integer values for actions
//	cout << AU << endl;
//	cout << AD << endl;
//	cout << AOU << endl;
//	cout << AOD << endl;
//	cout << AS << endl;
	//state->time_up[2] = 3;
	//state->time_up[4] = 3;
   	cout<<"0"<<endl;
   	cout<<flush; 
   	cerr<<"started1"; 
   	State present ;
   	// read from simulator  
    // update alight
    while(true)
    {   timeStamp++;
    	string s;  
	    getline(cin,s); 
	   	present.update(s);  
	   	//printS(present);
	   	nextStates.clear();
	   	present.elevator1.resetAlight();
	   	present.elevator2.resetAlight();
	   findNextStates2(1,1,present,1);   //TODO take care of the state with probability 1
	   Action nextAct= findAction() ;

	   changes(present,nextAct.first,nextAct.second) ;
	   
	   cout<<decodes[nextAct.first]<<"1 "<<decodes[nextAct.second]<<"2 "<<endl; 
		//cout<<nextStates.size();   		
		///printS(nextStates[5].first);
		//cout<<"prob"<<nextStates[5].second<<endl;
    }
    //
//	cout << state->toString() << endl;

//	int* distance1 = state->getDistanceArr(state->elevator1, AU);
//	int* distance2_down = state->getDistanceArr(state->elevator2, AD);
//	int* distance2_open = state->getDistanceArr(state->elevator2, AOU);
//	cout << "DOWN:" << endl;
//	state->printDistances(state->getMinDistanceArr(distance1, distance2_down));
//	cout << (state->getMinCost(distance1, distance2_down) /*+ state->insideCost(state->elevator1, distance1) + state->insideCost(state->elevator2, distance2_down)*/) << endl;
//	cout << "OPEN UP:" << endl;
//	state->printDistances(state->getMinDistanceArr(distance1, distance2_open));
//	cout << (state->getMinCost(distance1, distance2_open) /*+ state->insideCost(state->elevator1, distance1) + state->insideCost(state->elevator2, distance2_open)(*/) << endl;
//
//	Action action = state->getPolicyAction();
//	cout << "ACTIONS:" << endl;
//	cout << action.first << " " << action.second << endl << endl;
//
//	cout << state->toString() << endl;

	//double cost = state->runSimulation(1);

}