#include <iostream>
#include <limits>
using namespace std;

const int k = 2;
const int MAXN = 1024;
const int Bsize = 2 * k;  

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

// recursive function for insertion
node* insertRecursive(node* root, int point[], int disc) {
    if(root == NULL) {
        node* temp = new node;
        for(int i = 0; i < k; i++)
            temp->point[i] = point[i];
        temp->disc = disc;
        temp->left = temp->right = NULL;
        return temp;
    }

    if(arePointsSame(root->point, point))
        return root;

    int cd = root->disc;

    if(point[cd] < root->point[cd])
        root->left = insertRecursive(root->left, point, (cd + 1) % k);
    else
        root->right = insertRecursive(root->right, point, (cd + 1) % k);

    return root;
}

// function to insert a point(k dimensional array) in the tree
node* insert(node* root, int point[]){
    return insertRecursive(root, point, 0);
}

// recursive function for search
bool searchRecursive(node* root, int point[]) {
    if(root == NULL) return false;
    if(arePointsSame(root->point, point)) return true;

    int cd = root->disc;

    if(point[cd] < root->point[cd])
        return searchRecursive(root->left, point);

    return searchRecursive(root->right, point);
}

bool search(node* root, int point[]){
    return searchRecursive(root, point);
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
void regionSearch(node *p, double b[], double RECDEF[]){
    if(p == nullptr) return;

    if(inRegion(p, RECDEF))
        found(p);

    int cd = p->disc;

    double leftB[2*k], rightB[2*k];
    copyBounds(b, leftB, 2*k);
    copyBounds(b, rightB, 2*k);

    leftB[2*cd + 1] = p->point[cd];
    rightB[2*cd] = p->point[cd];

    if(p->left && boundsIntersectRegion(leftB, RECDEF))
        regionSearch(p->left, leftB, RECDEF);

    if(p->right && boundsIntersectRegion(rightB, RECDEF))
        regionSearch(p->right, rightB, RECDEF);
}

void partialMatch(node* root, int q[]) {
    double RECDEF[2*k];

    for(int i = 0; i < k; i++){
        if(q[i] == -1){
            RECDEF[2*i]   = -1e9;
            RECDEF[2*i+1] =  1e9;
        } 
        else {
            RECDEF[2*i]   = q[i];
            RECDEF[2*i+1] = q[i];
        }
    }

    double b[2*k];
    for(int i = 0; i < 2*k; i++)
        b[i] = (i%2==0 ? -1e9 : 1e9);

    cout << "Partial match results:\n";
    regionSearch(root, b, RECDEF);
}

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
void nearestNeighborRecursive(node* root, int target[], node*& best, double& bestDist){
    if(root == nullptr) return;

    double d = distanceSquared(root->point, target);
    if(d < bestDist){
        bestDist = d;
        best = root;
    }

    int cd = root->disc;

    node* nearChild;
    node* farChild;

    if(target[cd] < root->point[cd]){
        nearChild = root->left;
        farChild  = root->right;
    } else {
        nearChild = root->right;
        farChild  = root->left;
    }

    nearestNeighborRecursive(nearChild, target, best, bestDist);

    double diff = double(target[cd]) - double(root->point[cd]);
    if(diff*diff < bestDist)
        nearestNeighborRecursive(farChild, target, best, bestDist);
}

node* nearestNeighbor(node* root, int target[]){
    if(root == nullptr) return nullptr;
    node* best = nullptr;
    double bestDist = numeric_limits<double>::infinity();
    nearestNeighborRecursive(root, target, best, bestDist);
    return best;
}

// Query all points within a radius (circle) from target[]
void radiusQueryRecursive(node* root, int target[], double radiusSquared){
    if(root == nullptr) return;

    double d = distanceSquared(root->point, target);
    if(d <= radiusSquared){
        cout << "(" << root->point[0] << ", " << root->point[1] << ")\n";
    }

    int cd = root->disc;

    if(target[cd] < root->point[cd]){
        radiusQueryRecursive(root->left, target, radiusSquared);

        double diff = target[cd] - root->point[cd];
        if(diff*diff <= radiusSquared)
            radiusQueryRecursive(root->right, target, radiusSquared);
    }
    else{
        radiusQueryRecursive(root->right, target, radiusSquared);

        double diff = target[cd] - root->point[cd];
        if(diff*diff <= radiusSquared)
            radiusQueryRecursive(root->left, target, radiusSquared);
    }
}

void radiusQuery(node* root, int target[], double radius){
    double radiusSquared = radius * radius;
    cout << "Points inside radius " << radius << " from (";
    for(int i = 0; i < k; i++){
        cout << target[i];
        if(i+1 < k) cout << ", ";
    }
    cout << "):\n";

    radiusQueryRecursive(root, target, radiusSquared);
}

node* findMin(node* root, int dim) {
    if (root == nullptr) return nullptr;

    int cd = root->disc;

    if (cd == dim) {
        if (root->left == nullptr) return root;
        return findMin(root->left, dim);
    }

    node* leftMin  = findMin(root->left,  dim);
    node* rightMin = findMin(root->right, dim);

    node* minNode = root;
    if (leftMin  && leftMin->point[dim]  < minNode->point[dim])  minNode = leftMin;
    if (rightMin && rightMin->point[dim] < minNode->point[dim])  minNode = rightMin;

    return minNode;
}

node* findMax(node* root, int dim) {
    if (root == nullptr) return nullptr;

    int cd = root->disc;

    if (cd == dim) {
        if (root->right == nullptr) return root;
        return findMax(root->right, dim);
    }

    node* leftMax  = findMax(root->left,  dim);
    node* rightMax = findMax(root->right, dim);

    node* maxNode = root;
    if (leftMax  && leftMax->point[dim]  > maxNode->point[dim])  maxNode = leftMax;
    if (rightMax && rightMax->point[dim] > maxNode->point[dim])  maxNode = rightMax;

    return maxNode;
}

node* deleteNode(node* root, int point[]) {
    if (root == nullptr) return nullptr;

    int cd = root->disc;

    if (arePointsSame(root->point, point)) {

        if (root->right) {
            node* rep = findMin(root->right, cd);
            for (int i = 0; i < k; ++i) root->point[i] = rep->point[i];
            root->right = deleteNode(root->right, rep->point);
        }
        else if (root->left) {
            node* rep = findMax(root->left, cd);
            for (int i = 0; i < k; ++i) root->point[i] = rep->point[i];
            root->left = deleteNode(root->left, rep->point);
        }
        else {
            delete root;
            return nullptr;
        }

        return root;
    }

    if (point[cd] < root->point[cd])
        root->left = deleteNode(root->left, point);
    else
        root->right = deleteNode(root->right, point);

    return root;
}

void swapPoints(int a[], int b[], int k) {
    for (int i = 0; i < k; i++) {
        int temp = a[i];
        a[i] = b[i];
        b[i] = temp;
    }
}

int findMedianIndex(int A[][10], int n, int dim, int k) {

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

node* optimiseTree(int A[][10], int n, int dim, int k) {
    if (n == 0)
        return nullptr;

    int medianIndex = findMedianIndex(A, n, dim, k);

    node* root = new node();
    for (int i = 0; i < k; i++)
        root->point[i] = A[medianIndex][i];
    root->disc = dim;

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

    int nextDim = (dim + 1) % k;
    root->left  = optimiseTree(AL, leftN,  nextDim, k);
    root->right = optimiseTree(AH, rightN, nextDim, k);

    return root;
}

void collectPoints(node* p, int A[][10], int &idx, int k) {
    if (!p) return;

    for (int i = 0; i < k; i++)
        A[idx][i] = p->point[i];
    idx++;

    collectPoints(p->left, A, idx, k);
    collectPoints(p->right, A, idx, k);
}

void optimise(node* &root) {
    if (!root) {
        cout << "Tree is empty, nothing to optimise.\n";
        return;
    }

    int A[1000][10];
    int n = 0;

    collectPoints(root, A, n, k);
    root = optimiseTree(A, n, 0, k);

    cout << "\nOptimised KD-tree successfully built.\n";
}

void printTreeLevelOrder(node* root) {
    if (root == nullptr){
        cout << "[Empty tree]\n";
        return;
    }

    node* Q[1000];
    int front = 0, rear = 0;

    Q[rear++] = root;

    while (front < rear) {
        int levelSize = rear - front;
        cout << "Level: ";

        for (int i = 0; i < levelSize; i++) {
            node* cur = Q[front++];

            cout << "(" 
                 << cur->point[0] << ", "
                 << cur->point[1] << ") [d=" 
                 << cur->disc << "]  ";

            if (cur->left)  Q[rear++] = cur->left;
            if (cur->right) Q[rear++] = cur->right;
        }
        cout << "\n";
    }

}

// driver
int main() {
    node* root = nullptr;

    // initial KD-tree points
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

    cout << "\nBuilding KD-tree with initial points...\n";
    for(int i = 0; i < nPoints; ++i)
        root = insert(root, points[i]);

    while(true) {
        cout << "Choose option:\n";
        cout << "1. Print tree using BFS\n";
        cout << "2. Insert point\n";
        cout << "3. Search point\n";
        cout << "4. Delete point\n";
        cout << "5. Region Search\n";
        cout << "6. Partial Match\n";
        cout << "7. Nearest Neighbor\n";
        cout << "8. Radius Query\n";
        cout << "9. Optimise tree\n";
        cout << "10. Exit\n";
        cout << "Choose option: ";

        int choice;
        cin >> choice;

        if(choice == 10) break;

        int pt[k];

        switch(choice) {
        
        case 1:     //PRINT TREE
            cout << "Nodes in the tree:\n";
            printTreeLevelOrder(root);
            break;

        case 2:     // INSERT
            cout << "Enter point to insert:\n";
            for(int d = 0; d < k; d++){
                cout << "  Coord " << d << ": ";
                cin >> pt[d];
            }
            root = insert(root, pt);
            cout << "Inserted.\n";
            printTreeLevelOrder(root);
            break;

        case 3:     // SEARCH
            cout << "Enter point to search:\n";
            for(int d = 0; d < k; d++){
                cout << "  Coord " << d << ": ";
                cin >> pt[d];
            }
            cout << (search(root, pt) ? "Found.\n" : "Not found.\n");
            break;

        case 4:     // DELETE
            cout << "Enter point to delete:\n";
            for(int d = 0; d < k; d++){
                cout << "  Coord " << d << ": ";
                cin >> pt[d];
            }
            root = deleteNode(root, pt);
            cout << "Deleted (if existed).\n";
            printTreeLevelOrder(root);
            break;

        case 5: {   // REGION SEARCH
            double RECDEF[2*k];
            cout << "Enter region bounds [low, high] for each dimension:\n";
            for(int i = 0; i < k; i++){
                cout << "Dimension " << i << " low = ";
                cin >> RECDEF[2*i];
                cout << "Dimension " << i << " high = ";
                cin >> RECDEF[2*i+1];
            }

            double b[2*k];
            for(int i = 0; i < 2*k; i++)
                b[i] = (i%2==0 ? -1e9 : 1e9);

            cout << "Points inside region:\n";
            regionSearch(root, b, RECDEF);
            break;
        }

        case 6: {   // PARTIAL MATCH
            int q[k];
            cout << "Partial match query:\n";
            for(int i = 0; i < k; i++){
                int flag;
                cout << "Specify dimension " << i << "? (1=yes, 0=wildcard): ";
                cin >> flag;

                if(flag){
                    cout << "  Enter value: ";
                    cin >> q[i];
                } else {
                    q[i] = -1;   // wildcard flag
                }
            }
            partialMatch(root, q);
            break;
        }

        case 7: {   // NEAREST NEIGHBOR
            cout << "Enter target point:\n";
            for(int i = 0; i < k; i++){
                cout << "  Coord " << i << ": ";
                cin >> pt[i];
            }
            node* nn = nearestNeighbor(root, pt);
            if(nn){
                cout << "Nearest neighbor = (";
                for(int i = 0; i < k; i++){
                    cout << nn->point[i];
                    if(i+1 < k) cout << ", ";
                }
                cout << ")\n";
            }
            break;
        }

        case 8: {   // RADIUS QUERY
            double R;
            cout << "Enter center point:\n";
            for(int i = 0; i < k; i++){
                cout << "  Coord " << i << ": ";
                cin >> pt[i];
            }
            cout << "Enter radius: ";
            cin >> R;

            cout << "Points inside radius:\n";
            radiusQuery(root, pt, R);
            break;
        }

        case 9: {
            optimise(root);
            printTreeLevelOrder(root);
            break;
        }

        default:
            cout << "Invalid choice.\n";
        }
    }

    return 0;
}
