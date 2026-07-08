#include <iostream>
#include <limits>
#include <cmath>
using namespace std;

/*To do
~~~ done ~~~ Structure to make a point in k dim space as a node
~~~ done ~~~ Insertion 
Randomly built kd tree
~~~ done ~~~ Searching (exact match queries)
~~~ done ~~~ Insertion query 
~~~ done ~~~ Exact match query (basically search)
~~~ done ~~~ Partial match query
~~~ done ~~~ Region query
    - in region
    - found
    - bounds intersect region
~~~done~~~ NN match query
~~~ done ~~~ Deletion 
~~~ done ~~~ Optimal tree
Application 
Area of further research~~ 
extra: ~~~done ~~~radius query, 
bounding box
*/

const int k = 2;

// structure of a node in k-d tree
struct node
{
    int point[k];
    int disc;
    node *left, *right;
};


// to check if two points are same
bool arePointsSame (int point1[], int point2[]){
    for(int i = 0; i < k; i++){
        if(point1[i] != point2[i]) return false;
    }
    return true;
}

// function to create nodes (return type -> node)
node* createNode(int arr[])
{
    node* temp = new node;

    for(int i = 0; i < k; i++)
        temp->point[i] = arr[i];

    temp -> left = temp -> right = NULL;
    return temp;
}

// recursive function for insertion
node* insertRecursive(node* root, int point[], unsigned depth){

    if(root == NULL) 
        return createNode(point);

    if(arePointsSame(root->point, point))
        return root;
    
    unsigned cd = depth % k; //cd is the abbreviation for current dimension

    if(point[cd] < (root -> point[cd]))
        root -> left = insertRecursive(root -> left, point, depth + 1);
    else
        root -> right = insertRecursive(root -> right, point, depth+1);

    return root;
}

// function to insert a point(k dimensional array) in the tree
node* insert(node* root, int point[]){
    return insertRecursive(root, point, 0);
}

// recursive function for search
bool searchRecursive(node* root, int point[], unsigned depth){

    if(root == NULL) return false;
    if(arePointsSame(root -> point, point)) return true;

    unsigned cd = depth % k;

    if (point[cd] < root -> point[cd])
        return searchRecursive(root -> left, point, depth + 1);

    return searchRecursive(root -> right, point, depth + 1);
}

// function to search for a point
bool search(node* root, int point[]){
    return searchRecursive(root, point, 0);
}

// helper function to print points inside region
void found(node *p){
    cout << "(";
    for(int i = 0; i < k; i++){
        cout << p->point[i];
        if(i+1 < k) cout << ", ";
    }
    cout << ")\n";
}

// checks if a specific point is inside region
bool inRegion(node *p, double RECDEF[]){
    for(int i = 0; i < k; i++){
        if(p->point[i] < RECDEF[2*i] || p->point[i] > RECDEF[2*i+1])
            return false;
    }
    return true;
}

// checks if bounds of the subtree intersect(partially lie inside) queried boundary
bool boundsIntersectRegion(double b[], double RECDEF[]){
    for(int i = 0; i < k; i++){
        if(b[2*i] > RECDEF[2*i+1] || b[2*i+1] < RECDEF[2*i])
            return false;
    }
    return true;
}

void copyBounds(double a[], double b[], int size){
    for(int i = 0; i < size; i++){
        b[i] = a[i];
    }    
}

// find all points inside given query region
void regionSearch(node *p, double b[], double RECDEF[], int depth){
    if(p == nullptr) return;

    if(inRegion(p, RECDEF))
        found(p);

    int cd = depth % k;
    double boundsLeft[2*k], boundsRight[2*k];
    copyBounds(b, boundsLeft, 2*k);
    copyBounds(b, boundsRight, 2*k);

    boundsLeft[2*cd + 1] = p -> point[cd];
    boundsRight[2*cd] = p -> point[cd];

    if(p->left && boundsIntersectRegion(boundsLeft, RECDEF))
        regionSearch(p->left, boundsLeft, RECDEF, depth + 1);

    if(p->right && boundsIntersectRegion(boundsRight, RECDEF))
        regionSearch(p->right, boundsRight, RECDEF, depth + 1);
}
// -------------------- INSERTION QUERY -------------------- //

// Wrapper for "insertion query":
//  - checks if point already exists
//  - if not, inserts it and prints what happened
node* insertionQuery(node* root, int point[]){
    if(search(root, point)){
        cout << "Point (";
        for(int i = 0; i < k; i++){
            cout << point[i];
            if(i+1 < k) cout << ", ";
        }
        cout << ") already exists in the tree.\n";
        return root;
    }

    root = insert(root, point);

    cout << "Inserted point: (";
    for(int i = 0; i < k; i++){
        cout << point[i];
        if(i+1 < k) cout << ", ";
    }
    cout << ")\n";

    return root;
}

// Check if a node matches the partial query
bool matchesPartial(node* p, int query[], bool specified[]){
    for(int i = 0; i < k; i++){
        if(specified[i] && p->point[i] != query[i])
            return false;
    }
    return true;
}
// Recursive partial match search in k-d tree
void partialMatchRecursive(node* root, int query[], bool specified[],
                           unsigned depth)
{
    if(root == nullptr) return;

    // If this node matches the pattern, print it
    if(matchesPartial(root, query, specified)){
        found(root);   // reuses your existing "found" printer
    }

    unsigned cd = depth % k;

    // If this coordinate is not specified (wildcard), we must search both sides
    if(!specified[cd]){
        partialMatchRecursive(root->left,  query, specified, depth + 1);
        partialMatchRecursive(root->right, query, specified, depth + 1);
    } else {
        // Coordinate is specified: prune like normal kd-tree search
        if(query[cd] < root->point[cd])
            partialMatchRecursive(root->left,  query, specified, depth + 1);
        else
            partialMatchRecursive(root->right, query, specified, depth + 1);
    }
}

// --- Nearest Neighbor (NN) search helpers ---

// squared Euclidean distance between two k-d points
double distanceSquared(const int a[], const int b[]){
    double dist = 0.0;
    for(int i = 0; i < k; i++){
        double diff = double(a[i]) - double(b[i]);
        dist += diff * diff;
    }
    return dist;
}

// Recursive NN search
void nearestNeighborRecursive(node* root, int target[], unsigned depth, node*& best, double& bestDist){
    if(root == nullptr) return;

    // compute distance to current node
    double d = distanceSquared(root->point, target);
    if(d < bestDist){
        bestDist = d;
        best = root;
    }

    unsigned cd = depth % k;

    // Decide which side to visit first
    node* nearChild = nullptr;
    node* farChild  = nullptr;
    if(target[cd] < root->point[cd]){
        nearChild = root->left;
        farChild = root->right;
    } else {
        nearChild = root->right;
        farChild = root->left;
    }

    // Recurse into the near side first
    nearestNeighborRecursive(nearChild, target, depth + 1, best, bestDist);

    // Check whether we need to visit the far side:
    // If distance to splitting plane could be less than bestDist
    double diff = double(target[cd]) - double(root->point[cd]);
    double planeDist = diff * diff; // squared distance to plane
    if(planeDist < bestDist){
        nearestNeighborRecursive(farChild, target, depth + 1, best, bestDist);
    }
}

// Wrapper that returns the nearest node pointer (or nullptr if tree empty)
node* nearestNeighbor(node* root, int target[]){
    if(root == nullptr) return nullptr;
    node* best = nullptr;
    double bestDist = numeric_limits<double>::infinity();
    nearestNeighborRecursive(root, target, 0, best, bestDist);
    return best;
}

// -------------------- RADIUS QUERY -------------------- //

// Query all points within a radius (circle) from target[]
void radiusQueryRecursive(node* root, int target[], double radiusSquared,
                          unsigned depth)
{
    if(root == nullptr) return;

    // Check if current node lies inside circle
    double d = distanceSquared(root->point, target);
    if(d <= radiusSquared){
        cout << "(";
        for(int i = 0; i < k; i++){
            cout << root->point[i];
            if(i+1 < k) cout << ", ";
        }
        cout << ")\n";
    }

    unsigned cd = depth % k;

    // Always explore the near side
    if(target[cd] < root->point[cd])
    {
        radiusQueryRecursive(root->left, target, radiusSquared, depth+1);

        // Check if far side intersects the circle
        double diff = target[cd] - root->point[cd];
        if(diff*diff <= radiusSquared)
            radiusQueryRecursive(root->right, target, radiusSquared, depth+1);
    }
    else
    {
        radiusQueryRecursive(root->right, target, radiusSquared, depth+1);

        double diff = target[cd] - root->point[cd];
        if(diff*diff <= radiusSquared)
            radiusQueryRecursive(root->left, target, radiusSquared, depth+1);
    }
}

// Wrapper
void radiusQuery(node* root, int target[], double radius){
    double radiusSquared = radius * radius;
    cout << "Points inside radius " << radius << " from (";
    for(int i = 0; i < k; i++){
        cout << target[i];
        if(i+1 < k) cout << ", ";
    }
    cout << "):\n";

    radiusQueryRecursive(root, target, radiusSquared, 0);
}
void partialMatchQuery(node* root, int query[], bool specified[]){
    cout << "Partial match results:\n";
    partialMatchRecursive(root, query, specified, 0);
}

node* findMin(node* root, int dim, int depth) {
    if(root == nullptr) return nullptr;

    int cd = depth % k;

    // If this node splits on the same dimension,
    // minimum must be in the LEFT subtree only
    if(cd == dim) {
        if(root->left == nullptr) return root;
        return findMin(root->left, dim, depth + 1);
    }

    // Otherwise, minimum can be in ANY subtree or root
    node* leftMin  = findMin(root->left, dim, depth + 1);
    node* rightMin = findMin(root->right, dim, depth + 1);

    node* minNode = root;

    if(leftMin && leftMin->point[dim] < minNode->point[dim])
        minNode = leftMin;
    if(rightMin && rightMin->point[dim] < minNode->point[dim])
        minNode = rightMin;

    return minNode;
}
node* findMax(node* root, int dim, int depth) {
    if(root == nullptr) return nullptr;

    int cd = depth % k;

    if(cd == dim) {
        if(root->right == nullptr) return root;
        return findMax(root->right, dim, depth + 1);
    }

    node* leftMax  = findMax(root->left, dim, depth + 1);
    node* rightMax = findMax(root->right, dim, depth + 1);

    node* maxNode = root;

    if(leftMax && leftMax->point[dim] > maxNode->point[dim])
        maxNode = leftMax;
    if(rightMax && rightMax->point[dim] > maxNode->point[dim])
        maxNode = rightMax;

    return maxNode;
}
node* deleteNode(node* root, int point[], int depth = 0) {
    if(root == nullptr)
        return nullptr;

    int cd = depth % k;

    // If this is the node to delete
    if(arePointsSame(root->point, point)) {

        // CASE 1 : Node has right subtree
        if(root->right) {
            node* rep = findMin(root->right, cd, depth + 1);
            
            // Copy replacement data into root
            for(int i = 0; i < k; i++)
                root->point[i] = rep->point[i];

            // Delete the replacement from right subtree
            root->right = deleteNode(root->right, rep->point, depth + 1);
        }

        // CASE 2 : No right subtree, use left subtree
        else if(root->left) {
            node* rep = findMax(root->left, cd, depth + 1);

            // Copy replacement data
            for(int i = 0; i < k; i++)
                root->point[i] = rep->point[i];

            // Delete replacement node from left subtree
            root->left = deleteNode(root->left, rep->point, depth + 1);
        }

        // CASE 3 : Leaf node
        else {
            delete root;
            return nullptr;
        }

        return root;
    }

    // Traverse the tree normally
    if(point[cd] < root->point[cd])
        root->left = deleteNode(root->left, point, depth + 1);
    else
        root->right = deleteNode(root->right, point, depth + 1);

    return root;
}

// ---- helper: swap two points in array ----
void swapPoints(int a[], int b[], int k) {
    for (int i = 0; i < k; i++) {
        int temp = a[i];
        a[i] = b[i];
        b[i] = temp;
    }
}

// ---- partial selection to find median by dim (no STL) ----
int findMedianIndex(int A[][10], int n, int dim, int k) {
    // simple O(n^2) selection sort until median index
    int mid = n / 2;
    for (int i = 0; i <= mid; i++) {
        int minIndex = i;
        for (int j = i+1; j < n; j++) {
            if (A[j][dim] < A[minIndex][dim])
                minIndex = j;
        }
        swapPoints(A[i], A[minIndex], k);
    }
    return mid;
}

// ---- OPTIMIZE: builds balanced kd-tree ----
node* OPTIMIZE(int A[][10], int n, int dim, int k) {
    if (n == 0)
        return nullptr;

    // 1. find median
    int medianIndex = findMedianIndex(A, n, dim, k);

    // 2. create node
    node* root = new node();
    for (int i = 0; i < k; i++)
        root->point[i] = A[medianIndex][i];
    root->disc = dim;

    // 3. split into left and right arrays
    int AL[1000][10], AH[1000][10];
    int leftN = 0, rightN = 0;

    for (int i = 0; i < n; i++) {
        if (i == medianIndex) continue;

        if (A[i][dim] < A[medianIndex][dim]) {
            for (int d = 0; d < k; d++)
                AL[leftN][d] = A[i][d];
            leftN++;
        } else {
            for (int d = 0; d < k; d++)
                AH[rightN][d] = A[i][d];
            rightN++;
        }
    }

    // 4. recursive build
    int nextDim = (dim + 1) % k;
    root->left  = OPTIMIZE(AL, leftN,  nextDim, k);
    root->right = OPTIMIZE(AH, rightN, nextDim, k);

    return root;
}

// driver
int main(){
    node* root = nullptr;
    int points[][k] = {
        {3, 6},
        {17, 15},
        {13, 15},
        {6, 12},
        {9, 1},
        {2, 7},
        {10, 19}
    };

    const int nPoints = sizeof(points) / sizeof(points[0]);
    for(int i = 0; i < nPoints; ++i)
        root = insert(root, points[i]);
  
    int insPoint[k];
    cout << "Enter point for insertion query:\n";
    for(int d = 0; d < k; d++){
        cout << "  Coordinate " << d << ": ";
        cin >> insPoint[d];
    }
    root = insertionQuery(root, insPoint);

    int pmQuery[k];
    bool specified[k];

    cout << "\nPartial match query:\n";
    for(int d = 0; d < k; d++){
        int flag;
        cout << "Is coordinate " << d
             << " specified? (1 = yes, 0 = wildcard): ";
        cin >> flag;

        specified[d] = (flag != 0);

        if(specified[d]){
            cout << "  Enter value for coordinate " << d << ": ";
            cin >> pmQuery[d];
        } else {
            pmQuery[d] = 0;
        }
    }

    partialMatchQuery(root, pmQuery, specified);

    int target[k] = {9, 2};
    node* nn = nearestNeighbor(root, target);
    int query2[k] = {9, 2};
double R = 10.0;

cout << "\nRadius Query:\n";
radiusQuery(root, query2, R);


    if(nn){
        cout << "Query point: (" << target[0] << ", " << target[1] << ")\n";
        cout << "Nearest neighbor: (";
        for(int i = 0; i < k; i++){
            cout << nn->point[i];
            if(i+1 < k) cout << ", ";
        }
        cout << ")\n";
        cout << "Squared distance = " << distanceSquared(nn->point, target) << "\n";
    } else {
        cout << "Tree is empty.\n";
    }


     int p1[] = {3, 6};
    int p2[] = {17, 15};
    int p3[] = {13, 15};
    int p4[] = {6, 12};

    root = insert(root, p1);
    root = insert(root, p2);
    root = insert(root, p3);
    root = insert(root, p4);

    int toDelete[] = {17, 15};
    root = deleteNode(root, toDelete);

    return 0;
}
