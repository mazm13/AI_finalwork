#ifndef NODE_H
#define NODE_H

#include <iostream>
#include "myJudge.h"
#include <cstring>
#include <cstdlib>
#include <cmath>
#include <ctime>

using namespace std;

const int TREE_SIZE = 5000000;
const int MAX_SIZE = 13;
const int TIME_LIMITATION = 3500;
const double C = 0.8;

//每个节点代表一个局面
struct Node {
	int x, y; //下的棋子的位置;
	int profit;
	int visit_time;
	//e.g. belongs = 2时，代表machine在(x, y)下了一个棋子，之后的局面便是this node
	int belongs;
	void NodeSet(int _x, int _y, int _belongs) {
		x = _x;
		y = _y;
		belongs = _belongs;
		profit = 0;
		visit_time = 0;
	}
	Node() {
		x = -1; 
		y = -1;
		profit = 0;
		visit_time = 0;
	}
};

Node NodeTree[TREE_SIZE];

struct State {
	int M, N;
	int start_time;
	int noX, noY;
	int boardState[MAX_SIZE][MAX_SIZE]; //备用棋盘
	int topState[MAX_SIZE];				//备用top

	int board[MAX_SIZE][MAX_SIZE]; //公共的棋盘；
	int top[MAX_SIZE];

	State(int _M, int _N, const int* _top, const int* _board, int _noX, int _noY) \
	: M(_M), N(_N), noX(_noX), noY(_noY) {
		start_time = clock();
		for(int i = 0; i < N; ++i)
			top[i] = topState[i] = _top[i];
		for(int i = 0; i < M; ++i)
			for(int j = 0; j < N; ++j)
				board[i][j] = boardState[i][j] = _board[i * N + j];

		for(int i = 0; i < TREE_SIZE; ++i)
			NodeTree[i].NodeSet(0, 0, 0);
		NodeTree[0].NodeSet(-1, -1, 1);	//根节点代表user在(-1, -1)处下了一个棋子，并用1表示
	}

	void add_piece(int x, int y, int belongs) {
		board[x][y] = belongs;
		top[y]--;
		if((y == noY)&&(top[y] - 1 == noY)) top[y]--;
	}


	bool isTerminal(int n) {
		int x = NodeTree[n].x;
		int y = NodeTree[n].y;
		int belongs = NodeTree[n].belongs;
		if((x == -1)||(y == -1)) return false;
		if(((belongs == 2)&&(machineWin(x, y, M, N, board))) ||
		   ((belongs == 1)&&(userWin(x, y, M, N, board))) || isTie(N, top))
			return true;
		return false;
	}

	bool canbeExpanded(int n) {
		for(int i = 1; i <= N; ++i)
			if(NodeTree[n * N + i].visit_time == 0) return true;
		return false;
	}

	int expand(int n) {
		int en = 0;
		int exp_node[MAX_SIZE] = {0};
		for(int i = 1; i <= N; ++i) 
			if(NodeTree[n * N + i].visit_time == 0) exp_node[en++] = i - 1;
		//en == 0，说明当前节点不可扩展
		if(en == 0) return -1;
		int index = rand() % en;
		int y = exp_node[index];
		int x = top[y] - 1;
		add_piece(x, y, 3 - NodeTree[n].belongs);
		int j = n * N + y + 1;
		if(j > TREE_SIZE) {
			cout << "j > TREE_SIZE!!!!!" << endl;
			return 0;
		}
		NodeTree[j].NodeSet(x, y, 3 - NodeTree[n].belongs);
		return j;
	}

	int bestChild(int n) {
		int bc;
		double bp = - 10000000;
		int sum_vt = 0;
		for(int i = 1; i < N; ++i) {
			int j = n * N + i;
			if(NodeTree[j].visit_time == 0) continue;
			sum_vt += NodeTree[j].visit_time;
			double cp = 1.0 * NodeTree[j].profit / NodeTree[j].visit_time + \
						C * sqrt(2.0 * log(NodeTree[n].visit_time) / NodeTree[j].visit_time);
			if(cp > bp) {
				bp = cp;
				bc = j;
			}
		}
		if(sum_vt == 0) cout << "in function bestChild(), sum_vt == 0" << endl;
		return bc;
	}

	//树策略，往下选择节点的时候，跟着改变棋局
	int treePolicy(int n) {
		while(!isTerminal(n)) {
			int ne = expand(n);
			if(ne != -1) return ne;
			else {
				n = bestChild(n);
				add_piece(NodeTree[n].x, NodeTree[n].y, NodeTree[n].belongs);
			}
		}
		return n;
	}

	int defaultPolicy(int n) {

		if(isTerminal(n)) {
			int ddp = 1;
			if(NodeTree[n].belongs == 1) ddp = -1;
			return ddp;
		}

		int belongs = NodeTree[n].belongs;
		int cnt = belongs - 1;
		int dp = (belongs == 2) ? 1 : -1;
		for( ; ; cnt++) {
			int label[MAX_SIZE] = { 0 };
			int nl = 0;
			for(int i = 0; i < N; ++i)
				if(top[i] > 0) label[nl++] = i;
			int y1 = label[rand() % nl];
			int x1 = top[y1] - 1;
			if(cnt % 2 == 1) {
				add_piece(x1, y1, 1);
				if(userWin(x1, y1, M, N, board)) return -dp;
			}
			else {
				add_piece(x1, y1, 2);
				if(machineWin(x1, y1, M, N, board)) return dp;
			}
			if(isTie(N, top)) return 0;
		}
	}

	void backUp(int n, int delta) {
		while(n > 0) {
			NodeTree[n].visit_time ++;
			NodeTree[n].profit += delta;
			delta = -delta;
			n = (int)((n - 1) / N);
		}
		NodeTree[0].visit_time++;
		NodeTree[0].profit += delta;

		for(int i = 0; i < M; ++i)
			for(int j = 0; j < N; ++j)
				board[i][j] = boardState[i][j];
		for(int i = 0; i < N; ++i)
			top[i] = topState[i];
	}

	Node UCTSearch() {
		while(clock() - start_time < TIME_LIMITATION) {
			int select_node = treePolicy(0);
			int delta = defaultPolicy(select_node);
			backUp(select_node, delta);
		}
		int bc0 = bestChild(0);
		cout << endl << NodeTree[bestChild(0)].y << endl;
		cout << bestChild(0) << endl;
		return NodeTree[bc0];
	}
};

#endif
