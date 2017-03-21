#include "AIShell.h"

AIShell::AIShell(int numCols, int numRows, bool gravityOn, int** gameState, Move lastMove)
{
    this->deadline=5000;
    this->numRows=numRows;
    this->numCols=numCols;
    this->gravityOn=gravityOn;
    this->gameState=gameState;
    this->lastMove=lastMove;
}

AIShell::~AIShell()
{
    //delete the gameState variable.
    for (int i =0; i<numCols; i++){
        delete [] gameState[i]; 
    }
    delete [] gameState;

}

// Deletes the board.  Used when
// deallocating tempGameStates.
void AIShell::deleteBoard(int** board)
{
    for (int i = 0; i < numCols; i++) {
        delete[] board[i];
    }
    delete[] board;
}

// Swaps the first element with either the biggest/smallest
// element or the element to be pruned.
void AIShell::sort(std::vector<std::pair<int, int> >& moves, int index)
{
    std::pair<int, int> temp = moves[index];
    moves[index] = moves[0];
    moves[0] = temp;
}

// Note the limit is set at 4800 so IDSearch can
// exit out through all levels of recursion before
// the 5 minute mark
bool AIShell::isUnderLimit(struct timeval& start) const
{
    double result;
    struct timeval current;
    gettimeofday(&current, NULL);
    
    result = (current.tv_sec - start.tv_sec) * 1000.0;
    result += (current.tv_usec - start.tv_usec) / 1000.0;
    
    return result < 4800;
}

// Returns whether the move is invalid or not.
// Used in IDSearch after the search goes
// over the limit.
bool AIShell::isInvalidMove(Move& m)
{
    return m.col == -1 || m.row == -1;
}

std::vector<std::pair<int, int> > AIShell::emptyBlock(int**& gameState)
{
    
    std::vector<std::pair<int, int> > movesList;
    std::pair<int, int> oneMove;

    for (int col = 0; col < numCols; col++)
    {
        for (int row = 0; row < numRows; row++)
        {
            if (gameState[col][row] == NO_PIECE)
            {
                //std::cout << "Adding piece: " << "(" << col << "," << row << ")" << std::endl;
                oneMove = std::make_pair(col, row);
                movesList.push_back(oneMove);
                if(gravityOn) {break;}
            }
        }
    }

    return movesList;
}

// Returns the location of all threats on gamestate
std::vector<std::pair<int, int> > AIShell::threats(int**& gamestate)
{
    std::vector<std::pair<int, int> > threats;
    std::pair<int, int> aMove;
    
    for(int col = 0; col < numCols; col++)
    {
        for(int row = 0; row < numRows; row++)
        {
            if(gamestate[col][row] == NO_PIECE)
            {
                if(checkRow(col, row, HUMAN_PIECE, gamestate))
                {
                    aMove = std::make_pair(col, row);
                    threats.push_back(aMove);
                }
                
                else if(checkCol(col, row, HUMAN_PIECE, gamestate))
                {
                    aMove = std::make_pair(col, row);
                    threats.push_back(aMove);
                }
                
                else if(checkDiagonalLeftTop(col, row, HUMAN_PIECE, gamestate))
                {
                    aMove = std::make_pair(col, row);
                    threats.push_back(aMove);
                }
                
                else if(checkDiagonalLeftBot(col, row, HUMAN_PIECE, gamestate))
                {
                    aMove = std::make_pair(col, row);
                    threats.push_back(aMove);
                }
                
                if(gravityOn) {break;}
            }
        }
    }
    
    return threats;
}

// Returns a winning move if there is one, else
// it returns Move(-1, -1)
Move AIShell::winMoves(int**& gamestate)
{
    for(int col = 0; col < numCols; col++)
    {
        for(int row = 0; row < numRows; row++)
        {
            if(gamestate[col][row] == NO_PIECE)
            {
                if(checkRow(col, row, AI_PIECE, gamestate))
                {
                    return Move(col, row);
                }
                
                else if(checkCol(col, row, AI_PIECE, gamestate))
                {
                    return Move(col, row);
                }
                
                else if(checkDiagonalLeftTop(col, row, AI_PIECE, gamestate))
                {
                    return Move(col, row);
                }
                
                else if(checkDiagonalLeftBot(col, row, AI_PIECE, gamestate))
                {
                    return Move(col, row);
                }
                
                if(gravityOn) {break;}
            }
        }
    }
    
    return Move(-1, -1);
}

// checkRow checks all possible (k-1) in the same row from the given
// col and row parameters.
// For example: when col=4 and row=0 (4:0 on the board) and k=4,
// then checkRow checks the sequences [1:0, 2:0, 3:0, 4:0], [2:0, 3:0, 4:0, 5:0],
// [3:0, 4:0, 5:0, 6:0] and [4:0, 5:0, 6:0, 7:0] to see if there is a possible
// threat on the board or there's a possible winning move on the board.
bool AIShell::checkRow(int& col, int& row, int checkPlayer, int**& gamestate)
{
    int dCol = col-k+1;
    int dStartCol = dCol;
    
    int counter = 0;
    int numPlayer = 0;
    
    for(int i = 0; i < k; i++)
    {
        if(dCol >= 0 && dCol < numCols)
        {
            while(counter < k && dCol < numCols)
            {
                if(gamestate[dCol][row] == checkPlayer)
                {
                    numPlayer++;
                    counter++;
                }
                
                else if(gamestate[dCol][row] == NO_PIECE) {counter++;}
                else {break;}
                
                dCol++;
            }
            
            if(numPlayer == k-1) {return true;}
        }
        
        counter = 0;
        numPlayer = 0;
        dCol = dStartCol+i+1;
    }
    
    return false;
}

bool AIShell::checkCol(int& col, int& row, int checkPlayer, int**& gamestate)
{
    int dRow = row-k+1;
    int dStartRow = dRow;
    
    int counter = 0;
    int numPlayer = 0;
    
    for(int i = 0; i < k; i++)
    {
        if(dRow >= 0 && dRow < numRows)
        {
            while(counter < k && dRow < numRows)
            {
                if(gamestate[col][dRow] == checkPlayer)
                {
                    numPlayer++;
                    counter++;
                }
                
                else if(gamestate[col][dRow] == NO_PIECE) {counter++;}
                else {break;}
                
                dRow++;
            }
            
            if(numPlayer == k-1) {return true;}
        }
        
        counter = 0;
        numPlayer = 0;
        dRow = dStartRow+i+1;
    }
    
    return false;
}

bool AIShell::checkDiagonalLeftTop(int& col, int& row, int checkPlayer, int**& gamestate)
{
    int dCol = col-k+1;
    int dStartCol = dCol;
    
    int dRow = row+k-1;
    int dStartRow = dRow;
    
    int counter = 0;
    int numPlayer = 0;
    
    for(int i = 0; i < k; i++)
    {
        if(dCol >= 0 && dCol < numCols && dRow >= 0 && dRow < numRows)
        {
            while(counter < k && dCol < numCols && dRow >= 0)
            {
                if(gamestate[dCol][dRow] == checkPlayer)
                {
                    numPlayer++;
                    counter++;
                }
                
                else if(gamestate[dCol][dRow] == NO_PIECE) {counter++;}
                else {break;}
                
                dRow--;
                dCol++;
            }
            
            if(numPlayer == k-1) {return true;}
        }
        
        counter = 0;
        numPlayer = 0;
        dRow = dStartRow-i-1;
        dCol = dStartCol+i+1;
    }
    
    return false;
}

bool AIShell::checkDiagonalLeftBot(int& col, int& row, int checkPlayer, int**& gamestate)
{
    int dCol = col-k+1;
    int dStartCol = dCol;
    
    int dRow = row-k+1;
    int dStartRow = dRow;
    
    int counter = 0;
    int numPlayer = 0;
    
    for(int i = 0; i < k; i++)
    {
        if(dCol >= 0 && dCol < numCols && dRow >= 0 && dRow < numRows)
        {
            while(counter < k && dCol < numCols && dRow < numRows)
            {
                if(gamestate[dCol][dRow] == checkPlayer)
                {
                    numPlayer++;
                    counter++;
                }
                
                else if(gamestate[dCol][dRow] == NO_PIECE) {counter++;}
                else {break;}
                
                dRow++;
                dCol++;
            }
            
            if(numPlayer == k-1) {return true;}
        }
        
        counter = 0;
        numPlayer = 0;
        dRow = dStartRow+i+1;
        dCol = dStartCol+i+1;
    }
    
    return false;
}

unsigned int AIShell::winWays(int**& state, int player, struct timeval& start)
{
    unsigned int result = 0;
    
    for(int i = 0; i < numCols; i++)
    {
        for(int j = 0; j < numRows; j++)
        {
            if(state[i][j] == player)
            {
                std::string dirList[8]= {"N", "S", "W", "E", "NW", "NE", "SW", "SE"};
                for (int p=0; p<=7; p++)
                {
                    if(canWin(state, i, j, 0, player, dirList[p], start) == k) {result++;}
                    if(!isUnderLimit(start)) {return result;}
                }
            }
            
            if(!isUnderLimit(start)) {return result;}
        }
        
        if(!isUnderLimit(start)) {return result;}
    }
    
    return result;
}


bool AIShell::inDanger(int**& state, int colNum, int rowNum, int total, int player, std::string direction, struct timeval& start)
{
    int dangerCapacity = k-1;
    int hasSpace =0;
    if (direction == "N"){
        int tmpCol = colNum;
        int tmpRow = rowNum-1;
        while (tmpRow >=0 && isUnderLimit(start)){
            if (state[tmpCol][tmpRow]==player){
                return false;
            }
            if (state[tmpCol][tmpRow]!=NO_PIECE && state[tmpCol][tmpRow]!=player){
                dangerCapacity--;
            }
            if (state[tmpCol][tmpRow]==NO_PIECE){
                hasSpace++;
            }
            if (dangerCapacity<=0 && hasSpace<=2){
                return true;
            }
            tmpRow--;
        }
    }
    
    else if (direction == "S"){
        int tmpCol = colNum;
        int tmpRow = rowNum+1;
        while (tmpRow <= numRows-1 && isUnderLimit(start)){
            if (state[tmpCol][tmpRow]==player){
                return false;
            }
            if (state[tmpCol][tmpRow]!=NO_PIECE && state[tmpCol][tmpRow]!=player){
                dangerCapacity--;
            }
            if (state[tmpCol][tmpRow]==NO_PIECE){
                hasSpace++;
            }
            if (dangerCapacity<=0 && hasSpace<=2){
                return true;
            }
            tmpRow++;
        }
    }
    
    else if (direction == "W"){
        int tmpCol = colNum-1;
        int tmpRow = rowNum;
        while (tmpCol >=0 && isUnderLimit(start)){
            if (state[tmpCol][tmpRow]==player){
                return false;
            }
            if (state[tmpCol][tmpRow]!=NO_PIECE && state[tmpCol][tmpRow]!=player){
                dangerCapacity--;
            }
            if (state[tmpCol][tmpRow]==NO_PIECE){
                hasSpace++;
            }
            if (dangerCapacity<=0 && hasSpace<=2){
                return true;
            }
            tmpCol--;
        }
    }
    
    else if (direction == "E"){
        int tmpCol = colNum+1;
        int tmpRow = rowNum;
        while (tmpCol <= numCols-1 && isUnderLimit(start)){
            if (state[tmpCol][tmpRow]==player){
                return false;
            }
            if (state[tmpCol][tmpRow]!=NO_PIECE && state[tmpCol][tmpRow]!=player){
                dangerCapacity--;
            }
            if (state[tmpCol][tmpRow]==NO_PIECE){
                hasSpace++;
            }
            if (dangerCapacity<=0 && hasSpace<=2){
                return true;
            }
            tmpCol++;
        }
    }
    
    else if (direction == "NW"){
        int tmpCol = colNum-1;
        int tmpRow = rowNum-1;
        while (tmpCol >=0 && tmpRow >=0 && isUnderLimit(start)){
            if (state[tmpCol][tmpRow]==player){
                return false;
            }
            if (state[tmpCol][tmpRow]!=NO_PIECE && state[tmpCol][tmpRow]!=player){
                dangerCapacity--;
            }
            if (state[tmpCol][tmpRow]==NO_PIECE){
                hasSpace++;
            }
            if (dangerCapacity<=0 && hasSpace<=2){
                return true;
            }
            tmpCol--;
            tmpRow--;
        }
    }
    
    else if (direction == "NE"){
        int tmpCol = colNum+1;
        int tmpRow = rowNum-1;
        while (tmpCol <= numCols-1 && tmpRow >=0 && isUnderLimit(start)){
            if (state[tmpCol][tmpRow]==player){
                return false;
            }
            if (state[tmpCol][tmpRow]!=NO_PIECE && state[tmpCol][tmpRow]!=player){
                dangerCapacity--;
            }
            if (state[tmpCol][tmpRow]==NO_PIECE){
                hasSpace++;
            }
            if (dangerCapacity<=0 && hasSpace<=2){
                return true;
            }
            tmpCol++;
            tmpRow--;
        }
    }
    
    else if (direction == "SW"){
        int tmpCol = colNum-1;
        int tmpRow = rowNum+1;
        while (tmpCol >=0 && tmpRow <= numRows-1 && isUnderLimit(start)){
            if (state[tmpCol][tmpRow]==player){
                return false;
            }
            if (state[tmpCol][tmpRow]!=NO_PIECE && state[tmpCol][tmpRow]!=player){
                dangerCapacity--;
            }
            if (state[tmpCol][tmpRow]==NO_PIECE){
                hasSpace++;
            }
            if (dangerCapacity<=0 && hasSpace<=2){
                return true;
            }
            tmpCol--;
            tmpRow++;
        }
    }
    
    else if (direction == "SE"){
        int tmpCol = colNum+1;
        int tmpRow = rowNum+1;
        while (tmpCol <= numCols-1 && tmpRow <= numRows-1 && isUnderLimit(start)){
            if (state[tmpCol][tmpRow]==player){
                return false;
            }
            if (state[tmpCol][tmpRow]!=NO_PIECE && state[tmpCol][tmpRow]!=player){
                dangerCapacity--;
            }
            if (state[tmpCol][tmpRow]==NO_PIECE){
                hasSpace++;
            }
            if (dangerCapacity<=0 && hasSpace<=2){
                return true;
            }
            tmpCol++;
            tmpRow++;
        }
    }
    
    return false;
}


unsigned int AIShell::canWin(int**& state, int colNum, int rowNum, int total, int player, std::string direction, struct timeval& start)
{
    unsigned int result = 0;

    if(direction == "N")
    {
        for(int i = 0; i < k; i++)
        {
            if(rowNum+i >= numRows)
            {
                return result;
            }
            
            else
            {
                if(state[colNum][rowNum+i] == player || state[colNum][rowNum+i] == NO_PIECE)
                {
                    result++;
                }
                
                else
                {
                    return result;
                }
            }
        }
    }
    
    else if(direction == "S")
    {
        for(int i = 0; i < k; i++)
        {
            if(rowNum-i < 0)
            {
                return result;
            }
            
            else
            {
                if(state[colNum][rowNum-i] == player || state[colNum][rowNum-i] == NO_PIECE)
                {
                    result++;
                }
                
                else
                {
                    return result;
                }
            }
        }
    }
    
    else if(direction == "E")
    {
        for(int i = 0; i < k; i++)
        {
            if(colNum+i >= numCols)
            {
                return result;
            }
            
            else
            {
                if(state[colNum+i][rowNum] == player || state[colNum+i][rowNum] == NO_PIECE)
                {
                    result++;
                }
                
                else
                {
                    return result;
                }
            }
        }
    }
    
    else if(direction == "W")
    {
        for(int i = 0; i < k; i++)
        {
            if(colNum-i < 0)
            {
                return result;
            }
            
            else
            {
                if(state[colNum-i][rowNum] == player || state[colNum-i][rowNum] == NO_PIECE)
                {
                    result++;
                }
                
                else
                {
                    return result;
                }
            }
        }
    }
    
    else if(direction == "NW")
    {
        for(int i = 0; i < k; i++)
        {
            if(colNum-i < 0 || rowNum+i >= numRows)
            {
                return result;
            }
            
            else
            {
                if(state[colNum-i][rowNum+i] == player || state[colNum-i][rowNum+i] == NO_PIECE)
                {
                    result++;
                }
                
                else
                {
                    return result;
                }
            }
        }

    }
    
    else if(direction == "SE")
    {
        for(int i = 0; i < k; i++)
        {
            if(colNum+i >= numCols || rowNum-i < 0)
            {
                return result;
            }
            
            else
            {
                if(state[colNum+i][rowNum-i] == player || state[colNum+i][rowNum-i] == NO_PIECE)
                {
                    result++;
                }
                
                else
                {
                    return result;
                }
            }
        }
    }
    
    else if(direction == "NE")
    {
        for(int i = 0; i < k; i++)
        {
            if(colNum+i >= numCols || rowNum+i >= numRows)
            {
                return result;
            }
            
            else
            {
                if(state[colNum+i][rowNum+i] == player || state[colNum+i][rowNum+i] == NO_PIECE)
                {
                    result++;
                }
                
                else
                {
                    return result;
                }
            }
        }
    }
    
    else if(direction == "SW")
    {
        for(int i = 0; i < k; i++)
        {
            if(colNum-i < 0 || rowNum-i < 0)
            {
                return result;
            }
            
            else
            {
                if(state[colNum-i][rowNum-i] == player || state[colNum-i][rowNum-i] == NO_PIECE)
                {
                    result++;
                }
                
                else
                {
                    return result;
                }
            }
        }
    }
    
    return result;
}

int AIShell::heuristicEval(int** state, int player, struct timeval& start)
{
    int weight = 0;
    
    for(int i = 0; i < numCols; i++)
    {
        for(int j = 0; j < numRows; j++)
        {
            if(state[i][j] == player)
            {
                std::string dirList[8]= {"N", "S", "W", "E", "NW", "NE", "SW", "SE"};
                for (int p=0; p<=7; p++){
                    
                    // For each threat/win moves, add 100.
                    if(inDanger(state, i, j, 0, player, dirList[p], start))
                    {
                        weight += 100;
                    }
                    
                    if(!isUnderLimit(start)) {return 0;}
                }
            }
            
            if(!isUnderLimit(start)) {return 0;}
        }
        
        if(!isUnderLimit(start)) {return 0;}
    }

    if(player == HUMAN_PIECE)
    {
        weight *= -1;
    }
    
    return winWays(state, AI_PIECE, start) - winWays(state, HUMAN_PIECE, start) + weight;
}

std::pair<Move, int> AIShell::IDSearch()
{
    struct timeval start;
    gettimeofday(&start, NULL);
    
    int depth = 0;
    std::pair<Move, int> bestMove, temp;
    
    Move winMove = winMoves(gameState);
    if(!isInvalidMove(winMove))
    {
        bestMove = std::make_pair(winMove, 100);
        return bestMove;
    }
   
    do
    {
        depth++;
        temp = IDSearchRecurse(depth, gameState, AI_PIECE, -1000, 1000, start);
        
        if(!isInvalidMove(temp.first)){
            bestMove = temp;
        }
    }
    while(isUnderLimit(start));
    
    return bestMove;
}

std::pair<Move, int> AIShell::IDSearchRecurse(int depth, int** state, int turn, int alpha, int beta, struct timeval& start)
{
    std::pair<Move, int> inValid = std::make_pair(Move(-1, -1), -1);
    
    //AI_PIECE bigger better (MAX); HU_PIECE smaller better (MIN).
    int thisAlpha = alpha;
    int thisBeta = beta;
    
    // Create a copy of the current gameState
    int** tempGameState = new int*[numCols];
    for (int col = 0; col < numCols; col++)
    {
        tempGameState[col] = new int[numRows];
        for(int row = 0; row < numRows; row++)
        {
            tempGameState[col][row] = state[col][row];
        }
    }
    
    // Makes sure time is still under the limit
    if(!isUnderLimit(start)) {return inValid;}
    
    // Base Case for MinMax
    if (depth == 0)  //At leaf node, actual move is not important since we just care the eval of this node
    {
        Move dummyMove;
        int hValue = heuristicEval(state, turn, start);
        
        std::pair<Move, int> dummyMoveEval = std::make_pair(dummyMove, hValue);
        return dummyMoveEval;
    }
    
    // Recursive Step for MinMax
    else
    {
        // Search for threats (k-1 in a row) on the board
        std::vector<int> evalNodes;
        std::vector<std::pair<int, int> > movesList = threats(tempGameState);
        
        // If there are no threats find all available moves
        if(movesList.size() == 0 && isUnderLimit(start))
        {
            movesList = emptyBlock(tempGameState);
        }
      
        // Makes sure time is still under the limit
        if(!isUnderLimit(start)) {return inValid;}
        
        for (int i = 0; i < movesList.size(); i++)
        {
            // If there's only 1 possible move.
            if(movesList.size() == 1)
            {
                Move bestMove(movesList[i].first, movesList[i].second);
                std::pair<Move, int> evalValue = std::make_pair(bestMove, 0);
                return evalValue;
            }
            
            // If there's more than 1 possible move.
            tempGameState[movesList[i].first][movesList[i].second] = turn;
            
            if (turn == AI_PIECE){
                int nextLevelEval = IDSearchRecurse(depth - 1, tempGameState, HUMAN_PIECE, thisAlpha, thisBeta, start).second;
                // Makes sure time is still under the limit
                if(!isUnderLimit(start)) {return inValid;}
                
                evalNodes.push_back(nextLevelEval);
                if (nextLevelEval >= thisBeta){
                    sort(movesList, i);
                    evalNodes[0] = nextLevelEval;
                    break;
                }
                
                if (nextLevelEval >= thisAlpha){
                    sort(movesList, i);
                    evalNodes[0] = nextLevelEval;
                    thisAlpha=nextLevelEval;
                }
                
                // Makes sure time is still under the limit
                if(!isUnderLimit(start)) {return inValid;}
                
            }
            
            else{
                int nextLevelEval = IDSearchRecurse(depth - 1, tempGameState, AI_PIECE, thisAlpha, thisBeta, start).second;
                // Makes sure time is still under the limit
                    if(!isUnderLimit(start)) {return inValid;}
                
                evalNodes.push_back(nextLevelEval);
                if (nextLevelEval <= thisAlpha){
                    sort(movesList, i);
                    evalNodes[0] = nextLevelEval;
                    break;
                }
                
                if (nextLevelEval <= thisBeta){
                    sort(movesList, i);
                    evalNodes[0] = nextLevelEval;
                    thisBeta = nextLevelEval;
                }
                // Makes sure time is still under the limit
                if(!isUnderLimit(start)) {return inValid;}
            }
            
            tempGameState[movesList[i].first][movesList[i].second] = 0;   //Set the block back
            // Makes sure time is still under the limit
            if(!isUnderLimit(start)) {return inValid;}
        }
        
        if (turn == AI_PIECE)
        {
                
            Move bestMove(movesList[0].first, movesList[0].second);
            std::pair<Move, int> evalValue = std::make_pair(bestMove, evalNodes[0]);
                
            deleteBoard(tempGameState);
                
            // Makes sure time is still under the limit
            if(!isUnderLimit(start)) {return inValid;}
                
            return evalValue;
        }
            
        else
        {
            Move bestMove(movesList[0].first, movesList[0].second);
            std::pair<Move, int> evalValue = std::make_pair(bestMove, evalNodes[0]);
                
            deleteBoard(tempGameState);
                
            // Makes sure time is still under the limit
            if(!isUnderLimit(start)) {return inValid;}
                
            return evalValue;
        }
    }
}

Move AIShell::makeMove(){
    //this part should be filled in by the student to implement the AI
    //Example of a move could be: Move move(1, 2); //this will make a move at col 1, row 2
    
    //std::pair<Move, int> bestMove = makeMoveMinMax(3, gameState, AI_PIECE);
    //std::pair<Move, int> bestMove = makeMoveAlphaBeta(3, gameState, AI_PIECE, -1000, 1000);
    std::pair<Move, int> bestMove = IDSearch();
    
    return bestMove.first;
    
}
