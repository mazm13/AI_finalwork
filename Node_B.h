#ifndef NODE_H_
#define NODE_H_

#include <iostream>
#include <ctime>
#include "Judge.h"
#include <cstdlib>

using namespace std;

double C = 0.8;
const int TimeLimatation = 2000;
const int TREE_SIZE = 1000000;
const int MAX_SIZE = 12;

struct chess {
	int **board;
	int top[12];
	int noX, noY;
	int M, N;
	chess(int _M, int _N, const int _top[], const int _board[][12], int _noX, int _noY) {
		M = _M;
		N = _N;
		noX = _noX;
		noY = _noY;
		for(int i = 0; i < N; ++i)
			top[i] = _top[i];
		board = new int*[M];
		for(int i = 0; i < M; ++i) {
			board[i] = new int[N];
			for(int j = 0; j < N; ++j)
				board[i][j] = _board[i][j];
		}
	}
	bool User_win(int x, int y) {
		return userWin(x, y, M, N, board);
	}
	bool Machine_win(int x, int y) {
		return machineWin(x, y, M, N, board);
	}
	bool draw() {
		for(int i = 0; i < N; ++i)
			if(top[i] != 0) return false;
		return true;
	}
	bool user_add_piece(int y) {
		int _x = top[y] - 1;
		int _y = y;
		board[_x][_y] = 1;

		top[y]--;
		if((_x != 0)&&(y == noY)&&(_x - 1 == noX)) top[y]--;

		if(User_win(_x, _y)) return true;
		else return false;
	}
	bool machine_add_piece(int y) {
		int _x = top[y] - 1;
		int _y = y;
		board[_x][_y] = 2;

		top[y]--;
		if((_x != 0)&&(y == noY)&&(_x - 1 == noX)) top[y]--;

		if(Machine_win(_x, _y)) return true;
		else return false;
	}
	~chess() {
		for(int i = 0; i < M; ++i)
			delete [] board[i];
		delete [] board;
	}
};

struct Node {
	int M, N;
	int board[12][12];
	int top[12];
	int x, y;
	int profit, visit_time;
	int belongs;
	int expandNum;
	int expNode[12];
	bool visited;

	Node(int _M, int _N, const int _board[][12], const int _top[], int _belongs, int _x, int _y) {
		belongs = _belongs;
		x = _x;
		y = _y;
		M = _M;
		N = _N;
		for(int i = 0; i < N; ++i)
			top[i] = _top[i];
		for(int i = 0; i < M; ++i)
			for(int j = 0; j < N; ++j)
				board[i][j] = _board[i][j];
		expandNum = 0;
		for(int i = 0; i < N; ++i)
			if(top[i] > 0) expNode[expandNum++] = i;
		profit = 0;
		visit_time = 0;
		visited = true;
	}

	Node() {
		visited = false;
	}
};

Node NodeTree[TREE_SIZE];

struct State
{
	int M; 
	int N;
	int board[12][12];
	int top[12];
	int noX, noY;
	
	int start_time;

	State(int _M, int _N, const int* _top, const int* _board, int _noX, int _noY) \
	: M(_M), N(_N), noX(_noX), noY(_noY) {
		
		start_time = 0;
		for(int i = 0; i < N; ++i)
			top[i] = _top[i];
		for(int i = 0; i < M; ++i)
			for(int j = 0; j < N; ++j)
				board[i][j] = _board[i * N + j];
		NodeTree[0] = Node(M, N, board, top, 2, -1, -1);
		
	}

	bool isTerminal(int n) {
		Node node = NodeTree[n];

		if((node.x == -1)||(node.y == -1)) return false;

		bool isT = false;
		int **tempBoard;
		tempBoard = new int*[M];
		for(int i = 0; i < M; ++i) {
			tempBoard[i] = new int[N];
			for(int j = 0; j < N; ++j)
				tempBoard[i][j] = node.board[i][j];
		}
		int *ttop = new int[N];
		for(int i = 0; i < N; ++i)
			ttop[i] = top[i];
		
		if(((node.belongs == 2)&&(machineWin(node.x, node.y, M, N, tempBoard))) ||
		   ((node.belongs == 1)&&(userWin(node.x, node.y, M, N, tempBoard))) || isTie(N, ttop))
			isT = true;

		for(int i = 0; i < M; ++i)
			delete [] tempBoard[i];
		delete [] tempBoard;

		return isT;
	}

	bool canbeExpanded(int n) {
		return NodeTree[n].expandNum > 0;
	}

	int expand(int fn) {

		int _board[12][12];
		int _top[12];
		memcpy(_board, NodeTree[fn].board, sizeof(NodeTree[fn].board));
		memcpy(_top, NodeTree[fn].top, sizeof(NodeTree[fn].top));

		Node newNode = NodeTree[fn];
		int index = rand() % newNode.expandNum;
		int ny = newNode.expNode[index];
		int nx = newNode.top[ny] - 1;
		_board[nx][ny] = newNode.belongs;
		_top[ny] --;
		if((ny == noY)&&(nx - 1 == noX)) _top[ny]--;
		swap(NodeTree[fn].expNode[index], NodeTree[fn].expNode[--NodeTree[fn].expandNum]);

		Node node(M, N, _board, _top, 3 - newNode.belongs, nx, ny);
		NodeTree[fn * N + index + 1] = node;
		return (fn * N + index + 1);
	}

	int bestChild(int fn) {
		int bc;
		double bestProfit = -RAND_MAX;
		for(int i = 1; i <= N; ++i) {
			int j = fn * N + i;
			if(!NodeTree[j].visited) continue;
			double childProfit = 1.0 * NodeTree[j].profit / NodeTree[j].visit_time + \
								 C * sqrt(2.0 * log(NodeTree[fn].visit_time) / NodeTree[j].visit_time);
			if(childProfit > bestProfit){
				bc = j;
				bestProfit = childProfit;
			}
		}
		return bc;
	}

	//n 是选择的节点在树数组中的序数
	int treePolicy(int n) {
		while(!isTerminal(n)) {
			if(canbeExpanded(n))
				return expand(n);
			else
				n = bestChild(n);
		}
		return n;
	}

	int defaultPolicy(int n) {
		chess pc(M, N, NodeTree[n].top, NodeTree[n].board, noX, noY);
		int cnt = 1, dp = -1;
		if(NodeTree[n].belongs == 1) {
			cnt = 0;
			dp = 1;
		}
		for( ; ; cnt++) {
			int label[12] = {0};
			int nl = 0;
			for(int i = 0; i < N; ++i) {
				if(pc.top[i] > 0) label[nl++] = i;
			}
			int y1 = label[rand() % nl];
			if((cnt % 2 == 0)&&(pc.machine_add_piece(y1))) return dp;
			if((cnt % 2 == 1)&&(pc.user_add_piece(y1)))    return -dp;
			if(pc.draw()) return 0;
		}
	}

	void backUp(int n, int delta) {
		while(n > 0) {
			NodeTree[n].visit_time++;
			NodeTree[n].profit += delta;
			delta = -delta;
			n = (n - 1) / N;
		}
	}

	Node UCTSearch() {
		while(clock() - start_time < TimeLimatation) {
			int select_node = treePolicy(0);
			int delta = defaultPolicy(select_node);
			backUp(select_node, delta);
		}
		return NodeTree[bestChild(0)];
	}
};

#endif
