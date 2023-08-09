#include "gpbart.h"
#include <iomanip>
#include<cmath>
#include <random>
#include <RcppArmadillo.h>
using namespace std;

// =====================================
// Statistics Function
// =====================================

arma::mat rotate_operator(){

        arma::vec theta_vec_ = arma::linspace(0.0,M_PI,20);
        double theta_ = theta_vec_(arma::randi(arma::distr_param(0,(theta_vec_.size()-1))));

        // double theta_ = arma::randu(arma::distr_param(0.0,M_PI)); // Uniform grid
        double sin_ = std::sin(theta_);
        double cos_ = std::cos(theta_);

        arma::mat rot_mat_ = arma::mat(2,2);
        rot_mat_(0,0) = cos_;
        rot_mat_(0,1) = -sin_;
        rot_mat_(0,1) = sin_;
        rot_mat_(1,1) = cos_;

        return  rot_mat_;
}


//[[Rcpp::export]]
arma::mat inv_rcpp(arma::mat mat_ ,double tau) {
        return arma::inv_sympd(mat_+arma::eye(mat_.n_rows,mat_.n_rows)*(1/tau));
}

//[[Rcpp::export]]
arma::mat pinv_rcpp(arma::mat mat_ ,double tau) {
        return arma::pinv(mat_+arma::eye(mat_.n_rows,mat_.n_rows)*(1/tau));
}

// [[Rcpp::export]]
double gamma_pdf(double x, double a, double b) {

        double gamma_fun = tgamma(a);
        if(isinf(gamma_fun)){
                return 0.0;
        } else {
                return (pow(x, a-1) * exp(-x*b)*pow(b,a)) / ( gamma_fun);
        }


}

// [[Rcpp::export]]
double r_gamma_pdf(double x, double a, double b) {

        return R::dgamma(x,a,1/b,false);

}

// [[Rcpp::export]]
void print_mat_subset(arma::mat X) {
        int n_rows = X.n_rows;
        int n_cols = X.n_cols;

        // print the first 5 rows and 5 columns
        for (int i = 0; i < n_rows; i++) {
                if (i >= 5) break; // only print first 5 rows
                for (int j = 0; j < n_cols; j++) {
                        if (j >= 5) break; // only print first 5 columns
                        Rcpp::Rcout <<  std::setw(10) << X(i, j) << " ";
                }
                Rcpp::Rcout <<  std::endl;
        }
}


// Calculating the log-density of a MVN(0, Sigma)
//[[Rcpp::export]]
double log_dmvn(arma::vec& x, arma::mat& Sigma){

        arma::mat L = arma::chol(Sigma ,"lower"); // Remove diagonal later
        arma::vec D = L.diag();
        double p = Sigma.n_cols;

        arma::vec z(p);
        double out;
        double acc;

        for(int ip=0;ip<p;ip++){
                acc = 0.0;
                for(int ii = 0; ii < ip; ii++){
                        acc += z(ii)*L(ip,ii);
                }
                z(ip) = (x(ip)-acc)/D(ip);
        }
        out = (-0.5*sum(square(z))-( (p/2.0)*log(2.0*M_PI) +sum(log(D)) ));


        return out;

};

// //[[Rcpp::export]]
arma::mat sum_exclude_col(arma::mat mat, int exclude_int){

        // Setting the sum matrix
        arma::mat m(mat.n_rows,1);

        if(exclude_int==0){
                m = sum(mat.cols(1,mat.n_cols-1),1);
        } else if(exclude_int == (mat.n_cols-1)){
                m = sum(mat.cols(0,mat.n_cols-2),1);
        } else {
                m = arma::sum(mat.cols(0,exclude_int-1),1) + arma::sum(mat.cols(exclude_int+1,mat.n_cols-1),1);
        }

        return m;
}



// Implementation of the constructor outside of the struct
modelParam::modelParam(arma::mat x_train_,
                       arma::mat x_train_gp_,
                       arma::vec y_,
                       arma::mat x_test_,
                       arma::mat x_test_gp_,
                       int n_tree_,
                       int node_min_size_,
                       double alpha_,
                       double beta_,
                       double tau_mu_,
                       double tau_,
                       double a_tau_,
                       double d_tau_,
                       double nu_,
                       int n_mcmc_,
                       int n_burn_,
                       bool stump_bool_)
{
        // Assign the variables
        x_train = x_train_;
        x_train_gp = x_train_gp_;
        y = y_;
        x_test = x_test_;
        x_test_gp = x_test_gp_;
        n_tree = n_tree_;
        node_min_size = node_min_size_;
        alpha = alpha_;
        beta = beta_;
        tau_mu = tau_mu_;
        tau = tau_;
        a_tau = a_tau_;
        d_tau = d_tau_;
        nu = nu_;
        n_mcmc = n_mcmc_;
        n_burn = n_burn_;

        // Initialising distance and phi vector
        phi_grid = {0.1, 0.5, 1.0, 2.0, 3.0, 4.0, 50.0};

        phi_mat = arma::mat(n_tree_,x_train_gp.n_cols);


        // Random initialisation of phi_mat
        for(int i = 0; i < phi_mat.n_rows; i ++ ) {
                for(int j = 0; j < phi_mat.n_cols; j ++ ) {
                        phi_mat(i,j) = arma::randi(arma::distr_param(0,phi_grid.size()-1));
                }
        }


        // Grow acceptation ratio
        move_proposal = arma::vec(5,arma::fill::zeros);
        move_acceptance = arma::vec(5,arma::fill::zeros);


        stump_bool = stump_bool_; // Checking if only restrict the model to stumps
}


void compute_distance_matrices(modelParam& data,
                               arma::field<arma::cube>& pow_dist_train,
                               arma::field<arma::cube>& pow_dist_test,
                               arma::field<arma::cube>& pow_dist_train_test){

        // Rcpp::Rcout <<  "Grid size: " << data.phi_grid.size() << endl;

        for(int phi = 0; phi < data.phi_grid.size(); phi++){

                arma::cube pow_dist_train_(data.x_train_gp.n_rows,data.x_train_gp.n_rows,data.x_train_gp.n_cols,arma::fill::zeros);
                arma::cube pow_dist_test_(data.x_test_gp.n_rows,data.x_test_gp.n_rows,data.x_test_gp.n_cols,arma::fill::zeros);
                arma::cube pow_dist_train_test_(data.x_train_gp.n_rows,data.x_test_gp.n_rows,data.x_train_gp.n_cols,arma::fill::zeros);


                // Iterating over all
                for(int d = 0 ; d< data.x_train_gp.n_cols; d++ ){

                        // Iterating over the training observations
                        for(int  i = 0 ; i < data.x_train_gp.n_rows ; i ++) {

                                for(int i_plus = (i+1); i_plus < data.x_train_gp.n_rows ; i_plus ++ ){
                                        pow_dist_train_(i,i_plus,d) = pow_dist_train_(i_plus,i,d) = (data.x_train_gp(i,d)-data.x_train_gp(i_plus,d))*(data.x_train_gp(i,d)-data.x_train_gp(i_plus,d))/(2*data.phi_grid(phi)*data.phi_grid(phi)); // Remember that the power is the same regardless x_i or x_j
                                }
                                for(int j = 0 ; j < data.x_test_gp.n_rows ; j++) {
                                        pow_dist_train_test_(i,j,d) = ((data.x_train_gp(i,d)-data.x_test_gp(j,d))*(data.x_train_gp(i,d)-data.x_test_gp(j,d)))/(2*data.phi_grid(phi)*data.phi_grid(phi)); // Still need to fill the last row of pow_dist_train_test
                                }

                        }
                        // Iterating over the test observations
                        for(int k = 0; k< data.x_test_gp.n_rows ; k++) {
                                for(int  k_plus = k+1; k_plus < data.x_test_gp.n_rows; k_plus++){
                                        pow_dist_test_(k,k_plus,d) = pow_dist_test_(k_plus,k,d) = (data.x_test_gp(k,d)-data.x_test_gp(k_plus,d))*(data.x_test_gp(k,d)-data.x_test_gp(k_plus,d))/(2*data.phi_grid(phi)*data.phi_grid(phi));
                                }
                        }

                }

                // Iterating with the already evaluated phi
                pow_dist_train(phi,0,0) = pow_dist_train_;
                pow_dist_train_test(phi,0,0) = pow_dist_train_test_;
                pow_dist_test(phi,0,0) = pow_dist_test_;
        }
        return;

}


// Initialising a node
Node::Node(modelParam &data){
        isLeaf = true;
        isRoot = true;
        left = NULL;
        right = NULL;
        parent = NULL;
        train_index = arma::zeros<arma::vec>(data.x_train.n_rows);
        test_index = arma::zeros<arma::vec>(data.x_test.n_rows) ;
        train_index.fill(-1);
        test_index.fill(-1);

        // Initialising those vectors
        for(int i = 0; i< data.x_train.n_rows; i ++ ) {
                train_index[i] = -1;
        }
        for(int i = 0; i< data.x_test.n_rows; i ++ ) {
                test_index[i] = -1;
        }


        var_split = -1;
        var_split_rule = 0.0;
        lower = 0.0;
        upper = 1.0;
        mu = 0.0;
        n_leaf = 0.0;
        log_likelihood = 0.0;
        new_phi_log_likelihood = 0.0;
        depth = 0;




}

Node::~Node() {
        if(!isLeaf) {
                delete left;
                delete right;
        }
}

// Initializing a stump
void Node::Stump(modelParam& data){

        // Changing the left parent and right nodes;
        left = this;
        right = this;
        parent = this;
        // n_leaf  = data.x_train.n_rows;

        // Updating the training index with the current observations
        for(int i=0; i<data.x_train.n_rows;i++){
                train_index[i] = i;
        }

        // Updating the same for the test observations
        for(int i=0; i<data.x_test.n_rows;i++){
                test_index[i] = i;
        }

}

void Node::addingLeaves(modelParam& data){

     // Create the two new nodes
     left = new Node(data); // Creating a new vector object to the
     right = new Node(data);
     isLeaf = false;

     // Modifying the left node
     left -> isRoot = false;
     left -> isLeaf = true;
     left -> left = left;
     left -> right = left;
     left -> parent = this;
     left -> var_split = this->var_split; // This can break the code take care here
     left -> var_split_rule = this->var_split_rule; // This can break the code take care here
     left -> lower = 0.0;
     left -> upper = 1.0;
     left -> mu = 0.0;
     left -> log_likelihood = 0.0;
     left -> new_phi_log_likelihood = 0.0;
     left -> n_leaf = 0.0;
     left -> depth = this->depth+1;
     left -> train_index = arma::zeros<arma::vec>(data.x_train.n_rows);
     left -> test_index = arma::zeros<arma::vec>(data.x_test.n_rows);
     left -> train_index.fill(-1);
     left -> test_index.fill(-1);

     right -> isRoot = false;
     right -> isLeaf = true;
     right -> left = right; // Recall that you are saving the address of the right node.
     right -> right = right;
     right -> parent = this;
     right -> var_split = this->var_split; //This might break my code change here
     right -> var_split_rule = this->var_split_rule;
     right -> lower = 0.0;
     right -> upper = 1.0;
     right -> mu = 0.0;
     right -> log_likelihood = 0.0;
     right -> new_phi_log_likelihood = 0.0;
     right -> n_leaf = 0.0;
     right -> depth = this->depth+1;
     right -> train_index = arma::zeros<arma::vec>(data.x_train.n_rows);
     right -> test_index = arma::zeros<arma::vec>(data.x_test.n_rows);
     right -> train_index.fill(-1);
     right -> test_index.fill(-1);


     return;

}

// Creating boolean to check if the vector is left or right
bool Node::isLeft(){
        return (this == this->parent->left);
}

bool Node::isRight(){
        return (this == this->parent->right);
}

// Sample var
void Node::sampleSplitVar(modelParam &data){

          // Sampling one index from 0:(p-1)
          if(data.x_train.n_cols==1){
                  var_split = 0;
          } else {
                  var_split = arma::randi(arma::distr_param(0,(data.x_train.n_cols-1)));
          }

}


// This functions will get and update the current limits for this current variable
void Node::getLimits(){

        // Creating  a new pointer for the current node
        Node* x = this;
        // Already defined this -- no?
        lower = 0.0;
        upper = 1.0;
        // First we gonna check if the current node is a root or not
        bool tree_iter = x->isRoot ? false: true;
        while(tree_iter){
                bool is_left = x->isLeft(); // This gonna check if the current node is left or not
                x = x->parent; // Always getting the parent of the parent
                tree_iter = x->isRoot ? false : true; // To stop the while
                if(x->var_split == var_split){
                        tree_iter = false ; // This stop is necessary otherwise we would go up til the root, since we are always update there is no prob.
                        if(is_left){
                                upper = x->var_split_rule;
                                lower = x->lower;
                        } else {
                                upper = x->upper;
                                lower = x->var_split_rule;
                        }
                }
        }
}


void Node::displayCurrNode(){

                Rcpp::Rcout <<  "Node address: " << this << std::endl;
                Rcpp::Rcout <<  "Node parent: " << parent << std::endl;

                Rcpp::Rcout <<  "Cur Node is leaf: " << isLeaf << std::endl;
                Rcpp::Rcout <<  "Cur Node is root: " << isRoot << std::endl;
                Rcpp::Rcout <<  "Cur The split_var is: " << var_split << std::endl;
                Rcpp::Rcout <<  "Cur The split_var_rule is: " << var_split_rule << std::endl;

                return;
}


void Node::deletingLeaves(){

     // Should I create some warn to avoid memoery leak
     //something like it will only delete from a nog?
     // Deleting
     delete left; // This release the memory from the left point
     delete right; // This release the memory from the right point
     left = this;  // The new pointer for the left become the node itself
     right = this; // The new pointer for the right become the node itself
     isLeaf = true;

     return;

}
// Getting the leaves (this is the function that gonna do the recursion the
//                      function below is the one that gonna initialise it)
void get_leaves(Node* x,  std::vector<Node*> &leaves_vec) {

        if(x->isLeaf){
                leaves_vec.push_back(x);
        } else {
                get_leaves(x->left, leaves_vec);
                get_leaves(x->right,leaves_vec);
        }

        return;

}



// Initialising a vector of nodes in a standard way
std::vector<Node*> leaves(Node* x) {
        std::vector<Node*> leaves_init(0); // Initialising a vector of a vector of pointers of nodes of size zero
        get_leaves(x,leaves_init);
        return(leaves_init);
}

// Sweeping the trees looking for nogs
void get_nogs(std::vector<Node*>& nogs, Node* node){
        if(!node->isLeaf){
                bool bool_left_is_leaf = node->left->isLeaf;
                bool bool_right_is_leaf = node->right->isLeaf;

                // Checking if the current one is a NOGs
                if(bool_left_is_leaf && bool_right_is_leaf){
                        nogs.push_back(node);
                } else { // Keep looking for other NOGs
                        get_nogs(nogs, node->left);
                        get_nogs(nogs, node->right);
                }
        }
}

// Creating the vectors of nogs
std::vector<Node*> nogs(Node* tree){
        std::vector<Node*> nogs_init(0);
        get_nogs(nogs_init,tree);
        return nogs_init;
}



// Initializing the forest
Forest::Forest(modelParam& data){

        // Creatina vector of size of number of trees
        trees.resize(data.n_tree);
        for(int  i=0;i<data.n_tree;i++){
                // Creating the stump for each tree
                trees[i] = new Node(data);
                // Filling up each stump for each tree
                trees[i]->Stump(data);
        }
}

// Function to delete one tree
// Forest::~Forest(){
//         for(int  i=0;i<trees.size();i++){
//                 delete trees[i];
//         }
// }

// Selecting a random node
Node* sample_node(std::vector<Node*> leaves_){

        // Getting the number of leaves
        int n_leaves = leaves_.size();
        // return(leaves_[std::rand()%n_leaves]);
        if((n_leaves == 0) || (n_leaves==1) ){
                return leaves_[0];
        } else {
                return(leaves_[arma::randi(arma::distr_param(0,(n_leaves-1)))]);
        }

}
// Grow a tree for a given rule
void grow(Node* tree, modelParam &data, arma::vec &curr_res, int t,
          arma::field<arma::cube>& pow_dist_train,
          arma::field<arma::cube>& pow_dist_test,
          arma::field<arma::cube>& pow_dist_train_test){

        // Getting the number of terminal nodes
        std::vector<Node*> t_nodes = leaves(tree) ;
        std::vector<Node*> nog_nodes = nogs(tree);

        // Selecting one node to be sampled
        Node* g_node = sample_node(t_nodes);

        // Store all old quantities that will be used or not
        double old_lower = g_node->lower;
        double old_upper = g_node->upper;
        int old_var_split = g_node->var_split;
        double old_var_split_rule = g_node->var_split_rule;

        // Calculate current tree log likelihood
        double tree_log_like = 0;

        // Calculating the whole likelihood for the tree ( Do I really need to do this here? I do not think so)
        // (Yes, I do if I'm going to work with the loglikelihood over the trees to compute \phi_{i,j})
        // ----------------------------------------------
        // -- Gonna generalise the code to avoid this ---
        // ----------------------------------------------

        for(int i = 0; i < t_nodes.size(); i++){
                // Rcpp::Rcout <<  "Error gpNodeLogLike" << endl;
                t_nodes[i]->gpNodeLogLike(data, curr_res, t,
                                          pow_dist_train,
                                          pow_dist_test,
                                          pow_dist_train_test);
                tree_log_like = tree_log_like + t_nodes[i]->log_likelihood;
        }

        // Updating the grown_node only (NEED TO UPDATE ALL AGAIN IN THE)
        // g_node->gpNodeLogLike(data,curr_res,t);
        // tree_log_like = tree_log_like + g_node->log_likelihood;
        // Rcpp::Rcout <<  "LogLike Node ok Grow" << endl;

        // Adding the leaves
        g_node->addingLeaves(data);

        // Selecting the var
        g_node-> sampleSplitVar(data);
        // Updating the limits
        g_node->getLimits();


        // Selecting a rule
        g_node->var_split_rule = (g_node->upper-g_node->lower)*arma::randu(arma::distr_param(0.0,1.0))+g_node->lower;

        // Create an aux for the left and right index
        int train_left_counter = 0;
        int train_right_counter = 0;

        int test_left_counter = 0;
        int test_right_counter = 0;

        // Updating the left and the right nodes
        for(int i = 0;i<data.x_train.n_rows;i++){
                if(g_node -> train_index[i] == -1 ){
                        g_node->left->n_leaf = train_left_counter;
                        g_node->right->n_leaf = train_right_counter;
                        break;
                }
                if(data.x_train(g_node->train_index[i],g_node->var_split)<g_node->var_split_rule){
                        g_node->left->train_index[train_left_counter] = g_node->train_index[i];
                        train_left_counter++;
                } else {
                        g_node->right->train_index[train_right_counter] = g_node->train_index[i];
                        train_right_counter++;
                }

        }


        // Updating the left and right nodes for the
        for(int i = 0;i<data.x_test.n_rows; i++){
                if(g_node -> test_index[i] == -1){
                        g_node->left->n_leaf_test = test_left_counter;
                        g_node->right->n_leaf_test = test_right_counter;
                        break;
                }
                if(data.x_test(g_node->test_index[i],g_node->var_split)<g_node->var_split_rule){
                        g_node->left->test_index[test_left_counter] = g_node->test_index[i];
                        test_left_counter++;
                } else {
                        g_node->right->test_index[test_right_counter] = g_node->test_index[i];
                        test_right_counter++;
                }
        }

        // If is a root node
        if(g_node->isRoot){
                g_node->left->n_leaf = train_left_counter;
                g_node->right->n_leaf = train_right_counter;
                g_node->left->n_leaf_test = test_left_counter;
                g_node->right->n_leaf_test = test_right_counter;
        }

        // Avoiding nodes lower than the node_min
        if((g_node->left->n_leaf<data.node_min_size) || (g_node->right->n_leaf<data.node_min_size) ){

                // Rcpp::Rcout <<  " NODES" << endl;
                // Returning to the old values
                g_node->var_split = old_var_split;
                g_node->var_split_rule = old_var_split_rule;
                g_node->lower = old_lower;
                g_node->upper = old_upper;
                g_node->deletingLeaves();
                return;
        }


        // Updating the loglikelihood for those terminal nodes
        // Rcpp::Rcout <<  "Calculating likelihood of the new node on left" << endl;
        g_node->left->gpNodeLogLike(data, curr_res,t,
                                    pow_dist_train,
                                    pow_dist_test,
                                    pow_dist_train_test);
        // Rcpp::Rcout <<  "Calculating likelihood of the new node on right" << endl;
        g_node->right->gpNodeLogLike(data, curr_res,t,
                                     pow_dist_train,
                                     pow_dist_test,
                                     pow_dist_train_test);
        // Rcpp::Rcout <<  "NodeLogLike ok again" << endl;


        // Calculating the prior term for the grow
        double tree_prior = log(data.alpha*pow((1+g_node->depth),-data.beta)) +
                log(1-data.alpha*pow((1+g_node->depth+1),-data.beta)) + // Prior of left node being terminal
                log(1-data.alpha*pow((1+g_node->depth+1),-data.beta)) - // Prior of the right noide being terminal
                log(1-data.alpha*pow((1+g_node->depth),-data.beta)); // Old current node being terminal

        // Getting the transition probability
        double log_transition_prob = log((0.3)/(nog_nodes.size()+1)) - log(0.3/t_nodes.size()); // 0.3 and 0.3 are the prob of Prune and Grow, respectively

        // Calculating the loglikelihood for the new branches
        double new_tree_log_like = - g_node->log_likelihood + g_node->left->log_likelihood + g_node->right->log_likelihood ;

        // Calculating the acceptance ratio
        double acceptance = exp(new_tree_log_like  + log_transition_prob + tree_prior);

        if(data.stump_bool){
                acceptance = acceptance*(-1);
        }

        // Keeping the new tree or not
        if(arma::randu(arma::distr_param(0.0,1.0)) < acceptance){
                // Do nothing just keep the new tree
                // Rcpp::Rcout <<  " ACCEPTED" << endl;
                data.move_acceptance(0)++;
        } else {
                // Returning to the old values
                g_node->var_split = old_var_split;
                g_node->var_split_rule = old_var_split_rule;
                g_node->lower = old_lower;
                g_node->upper = old_upper;
                g_node->deletingLeaves();
        }

        return;

}



// Grow a tree for a given rule
void grow_rotation(Node* tree, modelParam &data, arma::vec &curr_res, int t,
                   arma::field<arma::cube>& pow_dist_train,
                   arma::field<arma::cube>& pow_dist_test,
                   arma::field<arma::cube>& pow_dist_train_test){

        // Getting the number of terminal nodes
        std::vector<Node*> t_nodes = leaves(tree) ;
        std::vector<Node*> nog_nodes = nogs(tree);

        // Selecting one node to be sampled
        Node* g_node = sample_node(t_nodes);

        // Store all old quantities that will be used or not
        double old_lower = g_node->lower;
        double old_upper = g_node->upper;
        int old_var_split = g_node->var_split;
        double old_var_split_rule = g_node->var_split_rule;

        // Calculate current tree log likelihood
        double tree_log_like = 0;


        // Calculating the whole likelihood for the tree ( Do I really need to do this here? I do not think so)
        // (Yes, I do if I'm going to work with the loglikelihood over the trees to compute \phi_{i,j})
        // ----------------------------------------------
        // -- Gonna generalize the code to avoid this ---
        // ----------------------------------------------

        // // Calculating the whole likelihood fo the tree
        for(int i = 0; i < t_nodes.size(); i++){
                // Rcpp::Rcout <<  "Error gpNodeLogLike" << endl;
                t_nodes[i]->gpNodeLogLike(data, curr_res, t,
                                          pow_dist_train,
                                          pow_dist_test,
                                          pow_dist_train_test); // Do I need to do this?
                tree_log_like = tree_log_like + t_nodes[i]->log_likelihood;
        }
        // Updating the grown_node only (NEED TO UPDATE ALL AGAIN IN THE)
        // g_node->gpNodeLogLike(data,curr_res,t);
        // tree_log_like = tree_log_like + g_node->log_likelihood;

        // Rcpp::Rcout <<  "LogLike Node ok Grow" << endl;

        // Adding the leaves
        g_node->addingLeaves(data);
        // Updating the limits
        g_node->getLimits();

        // Creating the rotation dummy aux
        arma::mat rotated_coord(data.x_train.n_rows,2);
        arma::mat rotated_coord_test(data.x_test.n_rows,2);

        // Sample var
        arma::vec candidates = arma::regspace(0,(data.x_train.n_cols-1));
        arma::vec sample = arma::shuffle(candidates);
        sample = sample.subvec(0,1);

        // Rotating coordinations
        rotated_coord.col(0) = data.x_train.col(sample(0));
        rotated_coord.col(1) = data.x_train.col(sample(1));

        rotated_coord_test.col(0) = data.x_test.col(sample(0));
        rotated_coord_test.col(1) = data.x_test.col(sample(1));

        arma::mat rotate_operator_ = rotate_operator();
        rotated_coord  = rotated_coord*rotate_operator_;
        rotated_coord_test = rotated_coord_test*rotate_operator_;

        // Selecting the var
        int selected_rot_var =  arma::randi(arma::distr_param(0,1));
        double rotated_var_split_rule = arma::randu(arma::distr_param(min(rotated_coord.col(selected_rot_var)),max(rotated_coord.col(selected_rot_var))));


        // Create an aux for the left and right index
        int train_left_counter = 0;
        int train_right_counter = 0;

        int test_left_counter = 0;
        int test_right_counter = 0;

        // Updating the left and the right nodes
        for(int i = 0;i<data.x_train.n_rows;i++){
                if(g_node -> train_index[i] == -1 ){
                        g_node->left->n_leaf = train_left_counter;
                        g_node->right->n_leaf = train_right_counter;
                        break;
                }
                if(rotated_coord(g_node->train_index[i],selected_rot_var)<rotated_var_split_rule){
                        g_node->left->train_index[train_left_counter] = g_node->train_index[i];
                        train_left_counter++;
                } else {
                        g_node->right->train_index[train_right_counter] = g_node->train_index[i];
                        train_right_counter++;
                }

        }


        // Updating the left and right nodes for the
        for(int i = 0;i<data.x_test.n_rows; i++){
                if(g_node -> test_index[i] == -1){
                        g_node->left->n_leaf_test = test_left_counter;
                        g_node->right->n_leaf_test = test_right_counter;
                        break;
                }
                if(rotated_coord_test(g_node->test_index[i],selected_rot_var)< rotated_var_split_rule){
                        g_node->left->test_index[test_left_counter] = g_node->test_index[i];
                        test_left_counter++;
                } else {
                        g_node->right->test_index[test_right_counter] = g_node->test_index[i];
                        test_right_counter++;
                }
        }

        // If is a root node
        if(g_node->isRoot){
                g_node->left->n_leaf = train_left_counter;
                g_node->right->n_leaf = train_right_counter;
                g_node->left->n_leaf_test = test_left_counter;
                g_node->right->n_leaf_test = test_right_counter;
        }

        // Avoiding nodes lower than the node_min
        if((g_node->left->n_leaf<data.node_min_size) || (g_node->right->n_leaf<data.node_min_size) ){

                // Rcpp::Rcout <<  " NODES" << endl;
                // Returning to the old values
                g_node->var_split = old_var_split;
                g_node->var_split_rule = old_var_split_rule;
                g_node->lower = old_lower;
                g_node->upper = old_upper;
                g_node->deletingLeaves();
                return;
        }


        // Updating the loglikelihood for those terminal nodes
        // Rcpp::Rcout <<  "Calculating likelihood of the new node on left" << endl;
        g_node->left->gpNodeLogLike(data, curr_res,t,
                                    pow_dist_train,
                                    pow_dist_test,
                                    pow_dist_train_test);
        // Rcpp::Rcout <<  "Calculating likelihood of the new node on right" << endl;
        g_node->right->gpNodeLogLike(data, curr_res,t,
                                     pow_dist_train,
                                     pow_dist_test,
                                     pow_dist_train_test);
        // Rcpp::Rcout <<  "NodeLogLike ok again" << endl;


        // Calculating the prior term for the grow
        double tree_prior = log(data.alpha*pow((1+g_node->depth),-data.beta)) +
                log(1-data.alpha*pow((1+g_node->depth+1),-data.beta)) + // Prior of left node being terminal
                log(1-data.alpha*pow((1+g_node->depth+1),-data.beta)) - // Prior of the right noide being terminal
                log(1-data.alpha*pow((1+g_node->depth),-data.beta)); // Old current node being terminal

        // Getting the transition probability
        double log_transition_prob = log((0.3)/(nog_nodes.size()+1)) - log(0.3/t_nodes.size()); // 0.3 and 0.3 are the prob of Prune and Grow, respectively

        // Calculating the loglikelihood for the new branches
        double new_tree_log_like = - g_node->log_likelihood + g_node->left->log_likelihood + g_node->right->log_likelihood ;

        // Calculating the acceptance ratio
        double acceptance = exp(new_tree_log_like  + log_transition_prob + tree_prior);

        if(data.stump_bool){
                acceptance = acceptance*(-1);
        }

        // Keeping the new tree or not
        if(arma::randu(arma::distr_param(0.0,1.0)) < acceptance){
                // Do nothing just keep the new tree
                // Rcpp::Rcout <<  " ACCEPTED" << endl;
                data.move_acceptance(1)++;
        } else {
                // Returning to the old values
                g_node->var_split = old_var_split;
                g_node->var_split_rule = old_var_split_rule;
                g_node->lower = old_lower;
                g_node->upper = old_upper;
                g_node->deletingLeaves();
        }

        return;

}


// Pruning a tree
void prune(Node* tree, modelParam&data, arma::vec &curr_res,  int t,
           arma::field<arma::cube>& pow_dist_train,
           arma::field<arma::cube>& pow_dist_test,
           arma::field<arma::cube>& pow_dist_train_test){


        // Getting the number of terminal nodes
        std::vector<Node*> t_nodes = leaves(tree);

        // Can't prune a root
        if(t_nodes.size()==1){
                // Rcpp::Rcout <<  "Nodes size " << t_nodes.size() <<endl;
                t_nodes[0]->gpNodeLogLike(data, curr_res, t,
                                          pow_dist_train,
                                          pow_dist_test,
                                          pow_dist_train_test);
                return;
        }

        std::vector<Node*> nog_nodes = nogs(tree);

        // Selecting one node to be sampled
        Node* p_node = sample_node(nog_nodes);


        // Calculate current tree log likelihood
        // double tree_log_like = 0;

        // Calculating the whole likelihood for the tree ( Do I really need to do this here? I do not think so)
        // (Yes, I do if I'm going to work with the loglikelihood over the trees to compute \phi_{i,j})
        // ----------------------------------------------
        // -- Gonna generalize the code to avoid this ---
        // ----------------------------------------------

        // Calculating the whole likelihood fo the tree
        for(int i = 0; i < t_nodes.size(); i++){
                t_nodes[i]->gpNodeLogLike(data, curr_res,t,
                                          pow_dist_train,
                                          pow_dist_test,
                                          pow_dist_train_test);
                // tree_log_like = tree_log_like + t_nodes[i]->log_likelihood;
        }

        // Rcpp::Rcout <<  "Error C1" << endl;
        // Updating the loglikelihood of the selected pruned node
        p_node->gpNodeLogLike(data, curr_res,t,
                              pow_dist_train,
                              pow_dist_test,
                              pow_dist_train_test);
        // p_node->left->gpNodeLogLike(data,curr_res,t);
        // p_node->right->gpNodeLogLike(data,curr_res,t);
        // Rcpp::Rcout <<  "Error C2" << endl;

        // Getting the loglikelihood of the new tree
        double new_tree_log_like =  p_node->log_likelihood - (p_node->left->log_likelihood + p_node->right->log_likelihood);

        // Calculating the transition loglikelihood
        double transition_loglike = log((0.3)/(t_nodes.size())) - log((0.3)/(nog_nodes.size()));

        // Calculating the prior term for the grow
        double tree_prior = log(1-data.alpha*pow((1+p_node->depth),-data.beta))-
                log(data.alpha*pow((1+p_node->depth),-data.beta)) -
                log(1-data.alpha*pow((1+p_node->depth+1),-data.beta)) - // Prior of left node being terminal
                log(1-data.alpha*pow((1+p_node->depth+1),-data.beta));  // Prior of the right noide being terminal
                 // Old current node being terminal


        // Calculating the acceptance
        double acceptance = exp(new_tree_log_like + transition_loglike + tree_prior);

        if(arma::randu(arma::distr_param(0.0,1.0))<acceptance){
                p_node->deletingLeaves();
                data.move_acceptance(2)++;
        } else {
                // p_node->left->gpNodeLogLike(data, curr_res);
                // p_node->right->gpNodeLogLike(data, curr_res);
        }

        return;
}


// // Creating the change verb
void change(Node* tree, modelParam &data, arma::vec &curr_res, int t,
            arma::field<arma::cube>& pow_dist_train,
            arma::field<arma::cube>& pow_dist_test,
            arma::field<arma::cube>& pow_dist_train_test){


        // Getting the number of terminal nodes
        std::vector<Node*> t_nodes = leaves(tree) ;
        std::vector<Node*> nog_nodes = nogs(tree);

        // Selecting one node to be sampled
        Node* c_node = sample_node(nog_nodes);

        // Calculate current tree log likelihood
        // double tree_log_like = 0; // Need to uncomment this too

        if(c_node->isRoot){
                // Rcpp::Rcout <<  " THAT NEVER HAPPENS" << endl;
               c_node-> n_leaf = data.x_train.n_rows;
               c_node-> n_leaf_test = data.x_test.n_rows;
        }

        // Calculating the whole likelihood for the tree ( Do I really need to do this here? I do not think so)
        // (Yes, I do if I'm going to work with the loglikelihood over the trees to compute \phi_{i,j})
        // ----------------------------------------------
        // -- Gonna generalize the code to avoid this ---
        // ----------------------------------------------
        // double tree_log_like  = 0.0;
        // Rcpp::Rcout <<  " Change error on terminal nodes" << endl;
        // Calculating the whole likelihood fo the tree
        for(int i = 0; i < t_nodes.size(); i++){
                // Rcpp::Rcout <<  "Loglike error " << ed
                t_nodes[i]->gpNodeLogLike(data, curr_res,t,
                                          pow_dist_train,
                                          pow_dist_test,
                                          pow_dist_train_test);
                // tree_log_like = tree_log_like + t_nodes[i]->log_likelihood;
        }

        // Updating the loglike of the nodes that gonna be changed only
        // c_node->left->gpNodeLogLike(data,curr_res,t);
        // c_node->right->gpNodeLogLike(data,curr_res,t);

        // Rcpp::Rcout <<  " Other kind of error" << endl;
        // If the current node has size zero there is no point of change its rule
        if(c_node->n_leaf==0) {
                return;
        }

        // Storing all the old loglikelihood from left
        double old_left_log_like = c_node->left->log_likelihood;
        arma::mat old_left_dist_train = c_node->left->dist_node_train;
        arma::vec old_left_leaf_res = c_node->left->leaf_res;



        arma::vec old_left_train_index = c_node->left->train_index;
        c_node->left->train_index.fill(-1); // Returning to the original
        int old_left_n_leaf = c_node->left->n_leaf;


        // Storing all of the old loglikelihood from right;
        double old_right_log_like = c_node->right->log_likelihood;
        arma::mat old_right_dist_train = c_node->right->dist_node_train;
        arma::vec old_right_leaf_res = c_node->right->leaf_res;


        arma::vec old_right_train_index = c_node->right->train_index;
        c_node->right->train_index.fill(-1);
        int old_right_n_leaf = c_node->right->n_leaf;



        // Storing test observations
        arma::vec old_left_test_index = c_node->left->test_index;
        arma::vec old_right_test_index = c_node->right->test_index;
        c_node->left->test_index.fill(-1);
        c_node->right->test_index.fill(-1);

        int old_left_n_leaf_test = c_node->left->n_leaf_test;
        int old_right_n_leaf_test = c_node->right->n_leaf_test;

        // ========
        // Storing the old ones
        int old_var_split = c_node->var_split;
        int old_var_split_rule = c_node->var_split_rule;
        int old_lower = c_node->lower;
        int old_upper = c_node->upper;

        // Selecting the var
        c_node-> sampleSplitVar(data);

        // Updating the limits
        c_node->getLimits();
        // Selecting a rule
        c_node -> var_split_rule = (c_node->upper-c_node->lower)*arma::randu(arma::distr_param(0.0,1.0))+c_node->lower;
        // c_node -> var_split_rule = 0.0;


        // ==========

        // Create an aux for the left and right index
        int train_left_counter = 0;
        int train_right_counter = 0;

        int test_left_counter = 0;
        int test_right_counter = 0;


        // Updating the left and the right nodes
        for(int i = 0;i<data.x_train.n_rows;i++){
                // Rcpp::Rcout <<  " Train indexeses " << c_node -> train_index[i] << endl ;
                if(c_node -> train_index[i] == -1){
                        c_node->left->n_leaf = train_left_counter;
                        c_node->right->n_leaf = train_right_counter;
                        break;
                }
                // Rcpp::Rcout <<  " Current train index " << c_node->train_index[i] << endl;

                if(data.x_train(c_node->train_index[i],c_node->var_split)<c_node->var_split_rule){
                        c_node->left->train_index[train_left_counter] = c_node->train_index[i];
                        train_left_counter++;
                } else {
                        c_node->right->train_index[train_right_counter] = c_node->train_index[i];
                        train_right_counter++;
                }
        }



        // Updating the left and the right nodes
        for(int i = 0;i<data.x_test.n_rows;i++){

                if(c_node -> test_index[i] == -1){
                        c_node->left->n_leaf_test = test_left_counter;
                        c_node->right->n_leaf_test = test_right_counter;
                        break;
                }

                if(data.x_test(c_node->test_index[i],c_node->var_split)<c_node->var_split_rule){
                        c_node->left->test_index[test_left_counter] = c_node->test_index[i];
                        test_left_counter++;
                } else {
                        c_node->right->test_index[test_right_counter] = c_node->test_index[i];
                        test_right_counter++;
                }
        }

        // If is a root node
        if(c_node->isRoot){
                c_node->left->n_leaf = train_left_counter;
                c_node->right->n_leaf = train_right_counter;
                c_node->left->n_leaf_test = test_left_counter;
                c_node->right->n_leaf_test = test_right_counter;
        }


        if((c_node->left->n_leaf<data.node_min_size) || (c_node->right->n_leaf)<data.node_min_size){

                // Returning to the previous values
                c_node->var_split = old_var_split;
                c_node->var_split_rule = old_var_split_rule;
                c_node->lower = old_lower;
                c_node->upper = old_upper;

                // Returning to the old ones
                c_node->left->dist_node_train = old_left_dist_train;
                c_node->left->leaf_res = old_left_leaf_res;

                c_node->left->n_leaf = old_left_n_leaf;
                c_node->left->n_leaf_test = old_left_n_leaf_test;
                c_node->left->log_likelihood = old_left_log_like;
                c_node->left->train_index = old_left_train_index;
                c_node->left->test_index = old_left_test_index;

                // Returning to the old ones
                c_node->right->dist_node_train = old_right_dist_train;
                c_node->right->leaf_res = old_right_leaf_res;

                c_node->right->n_leaf = old_right_n_leaf;
                c_node->right->n_leaf_test = old_right_n_leaf_test;
                c_node->right->log_likelihood = old_right_log_like;
                c_node->right->train_index = old_right_train_index;
                c_node->right->test_index = old_right_test_index;

                return;
        }

        // Updating the new left and right loglikelihoods
        c_node->left->gpNodeLogLike(data,curr_res,t,
                                    pow_dist_train,
                                    pow_dist_test,
                                    pow_dist_train_test);
        c_node->right->gpNodeLogLike(data,curr_res,t,
                                     pow_dist_train,
                                     pow_dist_test,
                                     pow_dist_train_test);

        // Calculating the acceptance
        double new_tree_log_like =  - old_left_log_like - old_right_log_like + c_node->left->log_likelihood + c_node->right->log_likelihood;

        double acceptance = exp(new_tree_log_like);

        if(arma::randu(arma::distr_param(0.0,1.0))<acceptance){
                // Keep all the trees
                data.move_acceptance(3)++;
        } else {

                // Returning to the previous values
                c_node->var_split = old_var_split;
                c_node->var_split_rule = old_var_split_rule;
                c_node->lower = old_lower;
                c_node->upper = old_upper;

                // Returning to the old ones
                c_node->left->dist_node_train = old_left_dist_train;
                c_node->left->leaf_res = old_left_leaf_res;


                c_node->left->n_leaf = old_left_n_leaf;
                c_node->left->n_leaf_test = old_left_n_leaf_test;
                c_node->left->log_likelihood = old_left_log_like;
                c_node->left->train_index = old_left_train_index;
                c_node->left->test_index = old_left_test_index;

                // Returning to the old ones
                c_node->right->dist_node_train = old_right_dist_train;
                c_node->right->leaf_res = old_right_leaf_res;

                c_node->right->n_leaf = old_right_n_leaf;
                c_node->right->n_leaf_test = old_right_n_leaf_test;
                c_node->right->log_likelihood = old_right_log_like;
                c_node->right->train_index = old_right_train_index;
                c_node->right->test_index = old_right_test_index;

        }

        return;
}


// // Creating the change verb
void change_rotation(Node* tree, modelParam &data, arma::vec &curr_res, int t,
                     arma::field<arma::cube>& pow_dist_train,
                     arma::field<arma::cube>& pow_dist_test,
                     arma::field<arma::cube>& pow_dist_train_test){


        // Getting the number of terminal nodes
        std::vector<Node*> t_nodes = leaves(tree) ;
        std::vector<Node*> nog_nodes = nogs(tree);

        // Selecting one node to be sampled
        Node* c_node = sample_node(nog_nodes);

        // Calculate current tree log likelihood
        // double tree_log_like = 0; // Need to uncomment this too

        if(c_node->isRoot){
                // Rcpp::Rcout <<  " THAT NEVER HAPPENS" << endl;
                c_node-> n_leaf = data.x_train.n_rows;
                c_node-> n_leaf_test = data.x_test.n_rows;
        }

        // Calculating the whole likelihood for the tree ( Do I really need to do this here? I do not think so)
        // (Yes, I do if I'm going to work with the loglikelihood over the trees to compute \phi_{i,j})
        // ----------------------------------------------
        // -- Gonna generalize the code to avoid this ---
        // ----------------------------------------------
        // double tree_log_like = 0.0;
        // // Rcpp::Rcout <<  " Change error on terminal nodes" << endl;
        // // Calculating the whole likelihood fo the tree
        for(int i = 0; i < t_nodes.size(); i++){
                // Rcpp::Rcout <<  "Loglike error " << ed
                t_nodes[i]->gpNodeLogLike(data, curr_res,t,
                                          pow_dist_train,
                                          pow_dist_test,
                                          pow_dist_train_test);
                // tree_log_like = tree_log_like + t_nodes[i]->log_likelihood;
        }


        // Updating the loglike of the nodes that gonna be changed only
        // c_node->left->gpNodeLogLike(data,curr_res,t);
        // c_node->right->gpNodeLogLike(data,curr_res,t);

        // Rcpp::Rcout <<  " Other kind of error" << endl;
        // If the current node has size zero there is no point of change its rule
        if(c_node->n_leaf==0) {
                return;
        }

        // Storing all the old loglikelihood from left
        double old_left_log_like = c_node->left->log_likelihood;
        arma::mat old_left_dist_train = c_node->left->dist_node_train;
        arma::vec old_left_leaf_res = c_node->left->leaf_res;


        arma::vec old_left_train_index = c_node->left->train_index;
        c_node->left->train_index.fill(-1); // Returning to the original
        int old_left_n_leaf = c_node->left->n_leaf;


        // Storing all of the old loglikelihood from right;
        double old_right_log_like = c_node->right->log_likelihood;
        arma::mat old_right_dist_train = c_node->right->dist_node_train;
        arma::vec old_right_leaf_res = c_node->right->leaf_res;


        arma::vec old_right_train_index = c_node->right->train_index;
        c_node->right->train_index.fill(-1);
        int old_right_n_leaf = c_node->right->n_leaf;


        // Storing test observations
        arma::vec old_left_test_index = c_node->left->test_index;
        arma::vec old_right_test_index = c_node->right->test_index;
        c_node->left->test_index.fill(-1);
        c_node->right->test_index.fill(-1);

        int old_left_n_leaf_test = c_node->left->n_leaf_test;
        int old_right_n_leaf_test = c_node->right->n_leaf_test;


        // Storing the old ones
        int old_var_split = c_node->var_split;
        int old_var_split_rule = c_node->var_split_rule;
        int old_lower = c_node->lower;
        int old_upper = c_node->upper;


        // Creating the rotation dummy aux
        arma::mat rotated_coord(data.x_train.n_rows,2);
        arma::mat rotated_coord_test(data.x_test.n_rows,2);

        // Sample var
        arma::vec candidates = arma::regspace(0,(data.x_train.n_cols-1));
        arma::vec sample = arma::shuffle(candidates);
        sample = sample.subvec(0,1);

        // Rotating coordinations
        rotated_coord.col(0) = data.x_train.col(sample(0));
        rotated_coord.col(1) = data.x_train.col(sample(1));

        rotated_coord_test.col(0) = data.x_test.col(sample(0));
        rotated_coord_test.col(1) = data.x_test.col(sample(1));

        arma::mat rotate_operator_ = rotate_operator();
        rotated_coord  = rotated_coord*rotate_operator_;
        rotated_coord_test = rotated_coord_test*rotate_operator_;

        // Selecting the var
        int selected_rot_var =  arma::randi(arma::distr_param(0,1));
        double rotated_var_split_rule = arma::randu(arma::distr_param(min(rotated_coord.col(selected_rot_var)),max(rotated_coord.col(selected_rot_var))));

        // Create an aux for the left and right index
        int train_left_counter = 0;
        int train_right_counter = 0;

        int test_left_counter = 0;
        int test_right_counter = 0;


        // Updating the left and the right nodes
        for(int i = 0;i<data.x_train.n_rows;i++){
                // Rcpp::Rcout <<  " Train indexeses " << c_node -> train_index[i] << endl ;
                if(c_node -> train_index[i] == -1){
                        c_node->left->n_leaf = train_left_counter;
                        c_node->right->n_leaf = train_right_counter;
                        break;
                }
                // Rcpp::Rcout <<  " Current train index " << c_node->train_index[i] << endl;

                if(rotated_coord(c_node->train_index[i],selected_rot_var)<rotated_var_split_rule){
                        c_node->left->train_index[train_left_counter] = c_node->train_index[i];
                        train_left_counter++;
                } else {
                        c_node->right->train_index[train_right_counter] = c_node->train_index[i];
                        train_right_counter++;
                }
        }



        // Updating the left and the right nodes
        for(int i = 0;i<data.x_test.n_rows;i++){

                if(c_node -> test_index[i] == -1){
                        c_node->left->n_leaf_test = test_left_counter;
                        c_node->right->n_leaf_test = test_right_counter;
                        break;
                }

                if(rotated_coord_test(c_node->test_index[i],selected_rot_var)<rotated_var_split_rule){
                        c_node->left->test_index[test_left_counter] = c_node->test_index[i];
                        test_left_counter++;
                } else {
                        c_node->right->test_index[test_right_counter] = c_node->test_index[i];
                        test_right_counter++;
                }
        }

        // If is a root node
        if(c_node->isRoot){
                c_node->left->n_leaf = train_left_counter;
                c_node->right->n_leaf = train_right_counter;
                c_node->left->n_leaf_test = test_left_counter;
                c_node->right->n_leaf_test = test_right_counter;
        }


        if((c_node->left->n_leaf<data.node_min_size) || (c_node->right->n_leaf)<data.node_min_size){

                // Returning to the previous values
                c_node->var_split = old_var_split;
                c_node->var_split_rule = old_var_split_rule;
                c_node->lower = old_lower;
                c_node->upper = old_upper;

                // Returning to the old ones
                c_node->left->dist_node_train = old_left_dist_train;
                c_node->left->leaf_res = old_left_leaf_res;

                c_node->left->n_leaf = old_left_n_leaf;
                c_node->left->n_leaf_test = old_left_n_leaf_test;
                c_node->left->log_likelihood = old_left_log_like;
                c_node->left->train_index = old_left_train_index;
                c_node->left->test_index = old_left_test_index;

                // Returning to the old ones
                c_node->right->dist_node_train = old_right_dist_train;
                c_node->right->leaf_res = old_right_leaf_res;

                c_node->right->n_leaf = old_right_n_leaf;
                c_node->right->n_leaf_test = old_right_n_leaf_test;
                c_node->right->log_likelihood = old_right_log_like;
                c_node->right->train_index = old_right_train_index;
                c_node->right->test_index = old_right_test_index;

                return;
        }

        // Updating the new left and right loglikelihoods
        c_node->left->gpNodeLogLike(data,curr_res,t,
                                    pow_dist_train,
                                    pow_dist_test,
                                    pow_dist_train_test);
        c_node->right->gpNodeLogLike(data,curr_res,t,
                                     pow_dist_train,
                                     pow_dist_test,
                                     pow_dist_train_test);

        // Calculating the acceptance
        double new_tree_log_like =  - old_left_log_like - old_right_log_like + c_node->left->log_likelihood + c_node->right->log_likelihood;

        double acceptance = exp(new_tree_log_like);

        if(arma::randu(arma::distr_param(0.0,1.0))<acceptance){
                // Keep all the trees
                data.move_acceptance(4)++;
        } else {

                // Returning to the previous values
                c_node->var_split = old_var_split;
                c_node->var_split_rule = old_var_split_rule;
                c_node->lower = old_lower;
                c_node->upper = old_upper;

                // Returning to the old ones
                c_node->left->dist_node_train = old_left_dist_train;
                c_node->left->leaf_res = old_left_leaf_res;


                c_node->left->n_leaf = old_left_n_leaf;
                c_node->left->n_leaf_test = old_left_n_leaf_test;
                c_node->left->log_likelihood = old_left_log_like;
                c_node->left->train_index = old_left_train_index;
                c_node->left->test_index = old_left_test_index;

                // Returning to the old ones
                c_node->right->dist_node_train = old_right_dist_train;
                c_node->right->leaf_res = old_right_leaf_res;

                c_node->right->n_leaf = old_right_n_leaf;
                c_node->right->n_leaf_test = old_right_n_leaf_test;
                c_node->right->log_likelihood = old_right_log_like;
                c_node->right->train_index = old_right_train_index;
                c_node->right->test_index = old_right_test_index;

        }

        return;
}



// Calculating the Loglilelihood of a node
void Node::gpNodeLogLike(modelParam& data, arma::vec &curr_res , int t,
                         arma::field<arma::cube>& pow_dist_train,
                         arma::field<arma::cube>& pow_dist_test,
                         arma::field<arma::cube>& pow_dist_train_test){

        // Getting number of leaves in case of a root
        if(isRoot){
                // Updating the r_sum
                n_leaf = data.x_train.n_rows;
                n_leaf_test = data.x_test.n_rows;
        }

        // Case of an empty node
        if(train_index[0]==-1){
        // if(n_leaf < 100){
                n_leaf = 0;
                log_likelihood = -2000000; // Absurd value avoid this case
                // Rcpp::Rcout <<  "OOOPS something happened" << endl;
                return;
        }


        // Creating the matrix elements to be filled
        arma::mat dist_node_train_(n_leaf,n_leaf);
        arma::vec leaf_res_(n_leaf);

        // Train elements
        for(int i = 0; i < n_leaf;i++){

                dist_node_train_(i,i) = 1.0;

                // CREATING THE TRAINING MATRIX COMMUNICATION
                for(int  k = (i+1); k < n_leaf; k++){
                        double total = 0.0;

                        for(int j = 0 ; j < data.x_train_gp.n_cols; j++){
                                total += pow_dist_train(data.phi_mat(t,j),0,0)(train_index[i],train_index[k],j);
                        }
                        // Rcpp::Rcout <<  "error here your fool" << endl;
                        dist_node_train_(k,i) = dist_node_train_(i,k) = exp(-total);
                }
                leaf_res_(i) = curr_res(train_index[i]);

        }

        // Sotring the new quantities for the node
        dist_node_train = dist_node_train_;
        leaf_res = leaf_res_;

        // Calculating the residuals covariance
        arma::mat res_cov(n_leaf,n_leaf);
        res_cov = (1/data.nu)*dist_node_train+arma::eye(n_leaf,n_leaf)*(1/data.tau)+(1/data.tau_mu);

        // If is smaller then the node size still need to update the quantities;
        if(n_leaf < data.node_min_size){
                log_likelihood = -2000000; // Absurd value avoid this case
                // Rcpp::Rcout <<  "OOOPS something happened" << endl;
                return;
        }

        // Getting the log-likelihood;
        log_likelihood = log_dmvn(leaf_res,res_cov);

        return;

}


// Calculating the Loglilelihood of a node
void Node::gpNodeLogLikeNewPhi(modelParam& data, arma::vec &curr_res,
                               arma::rowvec &new_phi_vec, int d, int t,
                               arma::field<arma::cube>& pow_dist_train,
                               arma::field<arma::cube>& pow_dist_test,
                               arma::field<arma::cube>& pow_dist_train_test
                               ){

        // Getting number of leaves in case of a root
        if(isRoot){
                // Updating the r_sum
                n_leaf = data.x_train.n_rows;
                n_leaf_test = data.x_test.n_rows;
        }


        // Case of an empty node
        if(train_index[0]==-1){
                // if(n_leaf < 100){
                n_leaf = 0;
                Rcpp::Rcout <<  "OOOPS something happened" << endl;

                return;
        }

        // Creating the matrix elements to be filled
        arma::mat dist_node_train_(n_leaf,n_leaf);
        // Train elements
        for(int i = 0; i < n_leaf;i++){
                dist_node_train_(i,i) = 1.0;
                // CREATING THE TRAINING MATRIX COMMUNICATION
                for(int  k = (i+1); k < n_leaf; k++){
                        double total = 0.0;
                        // Rcpp::Rcout <<  " New phi vec value" << new_phi_vec(0) << endl;
                        total = (-1)*log(dist_node_train(i,k)) - pow_dist_train(data.phi_mat(t,d),0,0)(train_index[i],train_index[k],d) + pow_dist_train(new_phi_vec(d),0,0)(train_index[i],train_index[k],d);
                        dist_node_train_(k,i) = dist_node_train_(i,k) = exp(-total);
                }

        }

        // Sotring the new quantities for the node
        new_phi_dist_node_train = dist_node_train_;

        // Calculating the residuals covariance
        arma::mat new_phi_res_cov(n_leaf,n_leaf);
        new_phi_res_cov = (1/data.nu)*new_phi_dist_node_train+arma::eye(n_leaf,n_leaf)*(1/data.tau)+(1/data.tau_mu);

        // If is smaller then the node size still need to update the quantities;
        if(n_leaf < data.node_min_size){
                new_phi_log_likelihood = -2000000; // Absurd value avoid this case
                return;
        }

        // Getting the log-likelihood;
        new_phi_log_likelihood = log_dmvn(leaf_res,new_phi_res_cov);

        return;

}


// UPDATING MU ( NOT NECESSARY)
// void updateMu(Node* tree, modelParam &data){
//
//         // Getting the terminal nodes
//         std::vector<Node*> t_nodes = leaves(tree);
//         // Iterating over the terminal nodes and updating the beta values
//         for(int i = 0; i < t_nodes.size();i++){
//
//                 if(t_nodes[i]->n_leaf==0 ){
//                         /// Skip nodes that doesn't have any observation within terminal node
//                         Rcpp::Rcout <<  " SKIPPED" << endl;
//                         continue;
//                 }
//
//
//                 // Calculating elements exclusive to each predictor
//                 // arma::mat aux_precision_inv = arma::inv(t_nodes[i]->B_t.slice(j)*t_nodes[i]->B.slice(j)+((data.n_tree*data.tau_b(j))/data.tau)*data.P); // Correction with the number of trees
//                 arma::mat psi_r_inv = arma::inv(arma::eye(t_nodes[i]->n_leaf,t_nodes[i]->n_leaf)*(1/data.tau) + (1/data.nu)*t_nodes[i]->dist_node_train);
//                 double aux_precision = arma::accu(psi_r_inv)+data.tau_mu;
//
//                 double mu_mean = arma::as_scalar(arma::ones(t_nodes[i]->n_leaf).t()*psi_r_inv*t_nodes[i]->leaf_res)/aux_precision;
//                 double mu_var = sqrt(1/aux_precision);
//
//                 // Updating mu
//                 t_nodes[i]->mu = arma::randn()*mu_var+mu_mean;
//         }
// }




// Get the prediction
// (MOST IMPORTANT AND COSTFUL FUNCTION FROM GP-BART)
void getPredictions(Node* tree,
                    modelParam &data,
                    arma::vec& current_prediction_train,
                    arma::vec& current_prediction_test,
                    int t,
                    bool sample,
                    arma::vec& y_hat_var,
                    arma::vec& y_hat_test_var,
                    arma::field<arma::cube>& pow_dist_train,
                    arma::field<arma::cube>& pow_dist_test,
                    arma::field<arma::cube>& pow_dist_train_test){

        // Getting the current prediction
        vector<Node*> t_nodes = leaves(tree);
        for(int i = 0; i<t_nodes.size();i++){

                // Skipping empty nodes
                if(t_nodes[i]->n_leaf==0){
                        Rcpp::Rcout <<  " THERE ARE EMPTY NODES" << endl;
                        continue;
                }

                // Calculating all the components to get the GP-samples
                arma::vec leaf_y_hat(t_nodes[i]->n_leaf,arma::fill::zeros);
                // Rcpp::Rcout <<  "Size of leaf test " << t_nodes[i]->n_leaf_test <<endl;
                arma::vec leaf_y_hat_test(t_nodes[i]->n_leaf_test, arma::fill::zeros);
                arma::mat Lambda = (1/data.nu)*t_nodes[i]->dist_node_train+(1/data.tau_mu);
                // arma::mat Lambda_inv = arma::inv(Lambda);
                // Rcpp::Rcout <<  "Entering in th matrix inversion" << endl;
                arma::mat Lambda_diag_inv;

                Lambda_diag_inv = arma::inv(Lambda + arma::eye(t_nodes[i]->n_leaf,t_nodes[i]->n_leaf)*(1/data.tau));


                // Getting the mean for the training
                arma::mat gp_train_mean = Lambda.t()*Lambda_diag_inv*t_nodes[i]->leaf_res;
                arma::mat gp_train_cov = Lambda - Lambda.t()*Lambda_diag_inv*Lambda;


                if(sample){
                        // Rcpp::Rcout <<  "SAMPLED TRAIN!" << endl;
                        leaf_y_hat = arma::mvnrnd(gp_train_mean,gp_train_cov+arma::eye(gp_train_cov.n_rows,gp_train_cov.n_cols)*1e-10);
                } else {
                        leaf_y_hat = gp_train_mean;
                }
                // Rcpp::Rcout <<  " Sample train ok" << endl;

                // For the training samples
                for(int j = 0; j<data.x_train.n_rows; j++){

                        if((t_nodes[i]->train_index[j])==-1.0){
                                break;
                        }
                        current_prediction_train[t_nodes[i]->train_index[j]] = leaf_y_hat[j];
                        y_hat_var[t_nodes[i]->train_index[j]] = gp_train_cov.diag().at(j)/pow(data.n_tree,2);

                }

                if(t_nodes[i]->n_leaf_test == 0 ){
                        continue;
                }


                // Creating the y_hat_test
                // This function is more complicated because it involves calculate the test and train new components;
                arma::mat dist_node_train_test_(t_nodes[i]->n_leaf,t_nodes[i]->n_leaf_test,arma::fill::ones);
                arma::mat dist_node_test_(t_nodes[i]->n_leaf_test,t_nodes[i]->n_leaf_test,arma::fill::ones);

                // ================================
                // Save this for the sampling new G
                // ================================
                // I should only calculate this if the tree is accepted actually;
                // CREATING THE TRAIN_TEST MATRIX
                for(int k = 0; k < t_nodes[i]->n_leaf_test; k++){
                        dist_node_test_(k,k) = 1.0;

                        // Rcpp::Rcout <<  "Error on k_star" << endl;
                        // Iterating with the training samples now
                        for(int k_train = 0; k_train< t_nodes[i]->n_leaf; k_train++){
                                double train_test_sum = 0.0;
                                for(int j_p=0;j_p < data.x_train_gp.n_cols;j_p++){
                                        // Rcpp::Rcout <<  "Analsying the phi from train-test  " << data.phi_vec(j_p) << endl;
                                        train_test_sum += pow_dist_train_test(data.phi_mat(t,j_p),0,0)(t_nodes[i]->train_index[k_train],t_nodes[i]->test_index[k],j_p);
                                }
                                dist_node_train_test_(k_train,k) = exp(-train_test_sum);
                        }

                        // Rcpp::Rcout <<  "Error on k_star_star" << endl;
                        for(int j_test = k+1; j_test < t_nodes[i]->n_leaf_test; j_test ++ ){
                                double test_sum = 0.0;
                                for(int j_p_test = 0; j_p_test < data.x_test_gp.n_cols; j_p_test++){ // Could fix this j_p_test loop
                                        // Rcpp::Rcout <<  "Analsying the phi from test  " << data.phi_vec(j_p_test) << endl;
                                        test_sum += pow_dist_test(data.phi_mat(t,j_p_test),0,0)(t_nodes[i]->test_index[k],t_nodes[i]->test_index[j_test],j_p_test);
                                }
                                // Storing the final result for those matrix;
                                dist_node_test_(j_test,k) = dist_node_test_(k,j_test) = exp(-test_sum);

                        }


                }


                // Getting the test samples;
                // Rcpp::Rcout <<  " Inverse nu " << (1/data.nu) << endl;
                arma::mat Lambda_star = (1/data.nu)*dist_node_train_test_+1/(data.tau_mu);
                arma::mat Lambda_star_star = (1/data.nu)*dist_node_test_+1/(data.tau_mu);

                arma::vec gp_mean_pred_test = Lambda_star.t()*Lambda_diag_inv*t_nodes[i]->leaf_res;
                arma::mat gp_cov_pred_test = Lambda_star_star - Lambda_star.t()*Lambda_diag_inv*Lambda_star;

                // arma::mat sample_test = arma::randn<arma::mat>(t_nodes[i]->n_leaf_test);

                // Rcpp::Rcout <<  "Error on the second cholesky decomposition" << endl;
                // leaf_y_hat_test = arma::chol(gp_cov_pred_test+arma::eye(gp_cov_pred_test.n_rows,gp_cov_pred_test.n_cols)*1e-8,"lower")*sample_test + gp_mean_pred_test;

                if(sample){
                        // Rcpp::Rcout <<  "SAMPLED!" << endl;
                        leaf_y_hat_test = arma::mvnrnd(gp_mean_pred_test,gp_cov_pred_test+arma::eye(gp_cov_pred_test.n_rows,gp_cov_pred_test.n_cols)*1e-10);
                } else {
                        leaf_y_hat_test = gp_mean_pred_test;

                }

                // Regarding the test samples
                for(int j = 0; j< data.x_test.n_rows;j++){

                        if(t_nodes[i]->test_index[j]==-1){
                                break;
                        }

                        current_prediction_test[t_nodes[i]->test_index[j]] = leaf_y_hat_test[j];
                        y_hat_test_var[t_nodes[i]->test_index[j]] = gp_cov_pred_test.diag().at(j)/pow(data.n_tree,2);

                }

        }
}

// Updating \phi
void updatePhi(Node* tree,
                    modelParam &data,
                    arma::vec& current_residuals,
                    int t,
                    arma::field<arma::cube>& pow_dist_train,
                    arma::field<arma::cube>& pow_dist_test,
                    arma::field<arma::cube>& pow_dist_train_test){

        // Getting the current prediction
        vector<Node*> t_nodes = leaves(tree);
        // Defining the grid of phis
        // arma::vec phi_grid = {0.1, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 5.5, 6.0, 6.5, 7.0, 7.5, 8.0, 8.5, 9.0, 9.5, 10.0, 10.0, 11.0, 12.0, 13.0, 14.0, 15.0, 16.0, 17.0, 18.0, 19.0, 20.0, 25.0, 30.0, 50.0, 75.0, 100.0, 125.0};
        // arma::vec phi_grid = {0.1, 0.5, 1.0, 2.5, 3.0, 4.0, 5.0, 10.0, 25.0, 50.0};



        // For iterating for each dimension p
        for(int  i = 0; i<data.x_train_gp.n_cols;i++){

                // Rcpp::Rcout <<  "Error on the phi sample" << endl;
                int sample_index = arma::randi(arma::distr_param(0,(data.phi_grid.size()-1))); //phi math now brings the index

                // Rcpp::Rcout <<  " Sampled phi is " << proposal_phi << endl;
                // Rcpp::Rcout <<  "Error storing the new phi vector" << endl;
                // Creating a copy of the current phi grid
                arma::rowvec new_phi_vec = data.phi_mat.row(t);
                new_phi_vec(i) = sample_index;
                // Rcpp::Rcout <<  "Size of the new phi vec is : " << new_phi_vec.size() << endl;
                // Rcpp::Rcout <<  "Size of the old phi vec is : " << data.phi_vec.size() << endl;

                // Skipping calculations in case of the same index is proposed;
                if(sample_index==data.phi_mat(t,i)){
                        continue;
                }

                // Storing loglikelihood for the ratios
                double old_phi_loglikelihood = 0.0;
                double new_phi_loglikelihood = 0.0;

                // Storing all quantities from the old tree
                for(int l = 0; l<t_nodes.size(); l++) {
                        // Rcpp::Rcout <<  "Error on the old phi loglikelihood" << endl;
                        // t_nodes[l]->gpNodeLogLike(data,current_residuals,t); // If I'm calculating for all terminal nodes that will not be necessary!
                        old_phi_loglikelihood += t_nodes[l]->log_likelihood;
                        // Rcpp::Rcout <<  "Error on the new function of phi loglikelihood" << endl;
                        t_nodes[l]->gpNodeLogLikeNewPhi(data,current_residuals,new_phi_vec,i,t,
                                                        pow_dist_train,
                                                        pow_dist_test,
                                                        pow_dist_train_test);
                        new_phi_loglikelihood += t_nodes[l]->new_phi_log_likelihood;

                }

                // Rcpp::Rcout <<  "The error is not here" << endl;
                // Verifying if it gonna accept or reject the new tree
                double prior_log = log(0.7*r_gamma_pdf(data.phi_grid(sample_index), 5000, 100) + 0.3*r_gamma_pdf(data.phi_grid(sample_index),3,(2.5))) - log(0.7*r_gamma_pdf(data.phi_grid(data.phi_mat(t,i)),5000,100)+0.3*r_gamma_pdf(data.phi_grid(data.phi_mat(t,i)),3,2.5));

                // Rcpp::Rcout <<  "Prior log value: " << prior_log << endl;

                if(isnan(prior_log)|| isnan(new_phi_loglikelihood) || isnan(old_phi_loglikelihood) ){
                        Rcpp::Rcout <<  " Something is WRONG (with phi)" << endl;
                }

                double acceptance = exp(new_phi_loglikelihood-old_phi_loglikelihood + prior_log);
                // Rcpp::Rcout <<  endl;
                // Rcpp::Rcout <<  "This corresponds to variable j = " << i << endl;
                // Rcpp::Rcout <<  "Current value is: " << data.phi_mat(i) << endl;
                // Rcpp::Rcout <<  "Sampled value was: " << new_phi_vec(i) << endl;
                // // Rcpp::Rcout <<  "Sampled value (TRUE): " << proposal_phi << endl;
                // Rcpp::Rcout <<  "The prior component has value: " << prior_log << endl;
                // Rcpp::Rcout <<  "New phi loglikelihood was: " << new_phi_loglikelihood << endl;
                // Rcpp::Rcout <<  "Old phi loglikelihood was: " << old_phi_loglikelihood << endl;
                // Rcpp::Rcout <<  endl;

                // Testing the parameters
                if(arma::randu(arma::distr_param(0.0,1.0)) < acceptance){
                        data.phi_mat(t,i) = sample_index; // Remember the line is the predictor and the column the tree
                        for(int l = 0; l<t_nodes.size(); l++) {
                                t_nodes[l]->dist_node_train = t_nodes[l]->new_phi_dist_node_train; // Is there a smart way of doing it?
                                t_nodes[l]->log_likelihood = t_nodes[l]->new_phi_log_likelihood;
                        }
                }
                // Rcpp::Rcout <<  "Printing the new phi vec:  " << data.phi_vec << endl;

        }

        // Ending the function
        return;
}


// Updating the tau parameter
void updateTau(arma::vec &y_hat,
               modelParam &data){

        // Getting the sum of residuals square
        double tau_res_sq_sum = dot((y_hat-data.y),(y_hat-data.y));

        data.tau = R::rgamma((0.5*data.y.size()+data.a_tau),1/(0.5*tau_res_sq_sum+data.d_tau));

        return;
}


// Creating the BART function
// [[Rcpp::export]]
Rcpp::List cppgpbart(arma::mat x_train,
                     arma::mat x_train_gp,
                          arma::vec y_train,
                          arma::mat x_test,
                          arma::mat x_test_gp,
                          int n_tree,
                          int node_min_size,
                          int n_mcmc,
                          int n_burn,
                          double tau, double mu,
                          double tau_mu,
                          double alpha, double beta,
                          double a_tau, double d_tau,
                          double nu,
                          bool stump_bool,
                          bool sample,
                          bool sample_phi,
                          bool verbose,
                          bool no_rotation_bool,
                          bool only_rotation_bool){

        // Posterior counter
        int curr = 0;


        // Rcpp::Rcout <<  " Error on model.param" << endl;
        // Creating the structu object
        modelParam data(x_train,
                        x_train_gp,
                        y_train,
                        x_test,
                        x_test_gp,
                        n_tree,
                        node_min_size,
                        alpha,
                        beta,
                        tau_mu,
                        tau,
                        a_tau,
                        d_tau,
                        nu,
                        n_mcmc,
                        n_burn,
                        stump_bool);


        arma::field<arma::cube> pow_dist_train(10,1,1);
        arma::field<arma::cube> pow_dist_test(10,1,1);
        arma::field<arma::cube> pow_dist_train_test(10,1,1);


        // Rcpp::Rcout <<  "error here" << endl;
        // Computing all elements needed to the distance matrix
        compute_distance_matrices(data,
                                  pow_dist_train,
                                  pow_dist_test,
                                  pow_dist_train_test);

        // Rcpp::Rcout <<  "error computing matrices" << endl;


        // Getting the n_post
        int n_post = n_mcmc - n_burn;

        // Creating a vector to store the proposal for each

        // Defining those elements
        arma::mat y_train_hat_post = arma::zeros<arma::mat>(data.x_train.n_rows,n_post);
        arma::mat y_test_hat_post = arma::zeros<arma::mat>(data.x_test.n_rows,n_post);
        arma::mat y_train_var_post = arma::zeros<arma::mat>(data.x_train.n_rows,n_post);
        arma::mat y_test_var_post = arma::zeros<arma::mat>(data.x_test.n_rows,n_post);


        arma::cube all_tree_post(y_train.size(),n_tree,n_post,arma::fill::zeros);
        arma::cube all_phi_post(n_tree,x_train_gp.n_cols,n_post,arma::fill::zeros);
        arma::vec tau_post = arma::zeros<arma::vec>(n_post);
        arma::vec tau_post_all = arma::vec(n_mcmc,arma::fill::zeros);

        // Defining other variables
        // arma::vec partial_pred = arma::zeros<arma::vec>(data.x_train.n_rows);
        arma::vec partial_pred = (data.y)/n_tree;
        arma::vec partial_residuals = arma::zeros<arma::vec>(data.x_train.n_rows);
        arma::mat tree_fits_store(data.x_train.n_rows,data.n_tree,arma::fill::zeros);
        for(int i = 0 ; i < data.n_tree ; i ++ ){
                tree_fits_store.col(i) = partial_pred;
        }
        arma::mat tree_fits_store_test(data.x_test.n_rows,data.n_tree,arma::fill::zeros);

        // Variance
        arma::mat tree_var_store(data.x_train.n_rows,data.n_tree,arma::fill::zeros);
        arma::mat tree_var_test_store(data.x_test.n_rows,data.n_tree,arma::fill::zeros);



        double verb;

        // Defining progress bars parameters
        const int width = 70;
        double pb = 0;


        // Rcpp::Rcout <<  " Error one " << endl;

        // Selecting the train
        Forest all_forest(data);

        for(int i = 0;i<data.n_mcmc;i++){

                // Initialising PB
                if(verbose){

                        Rcpp::Rcout <<  "[";
                        int k = 0;
                        // Evaluating progress bar
                        for(;k<=pb*width/data.n_mcmc;k++){
                                Rcpp::Rcout <<  "=";
                        }

                        for(; k < width;k++){
                                Rcpp::Rcout <<  " ";
                        }

                        Rcpp::Rcout <<  "] " << std::setprecision(5) << (pb/data.n_mcmc)*100 << "%\r";
                        Rcpp::Rcout.flush();


                }

                // Getting zeros
                arma::vec prediction_train_sum(data.x_train.n_rows,arma::fill::zeros);
                arma::vec prediction_test_sum(data.x_test.n_rows,arma::fill::zeros);
                arma::vec var_train_sum(data.x_train.n_rows,arma::fill::zeros);
                arma::vec var_test_sum(data.x_test.n_rows,arma::fill::zeros);

                for(int t = 0; t<data.n_tree;t++){


                        // Creating the auxliar prediction vector
                        arma::vec y_hat(data.y.n_rows,arma::fill::zeros);
                        arma::vec prediction_test(data.x_test.n_rows,arma::fill::zeros);
                        arma::vec y_hat_var(data.y.n_rows,arma::fill::zeros);
                        arma::vec y_hat_test_var(data.x_test.n_rows,arma::fill::zeros);



                        // Rcpp::Rcout <<  "Residuals error "<< endl;
                        // Updating the partial residuals
                        if(data.n_tree>1){
                                partial_residuals = data.y-sum_exclude_col(tree_fits_store,t);

                        } else {
                                partial_residuals = data.y;
                        }

                        // Iterating over all trees
                        verb = arma::randu(arma::distr_param(0.0,1.0));
                        if(all_forest.trees[t]->isLeaf & all_forest.trees[t]->isRoot){
                                verb = arma::randu(arma::distr_param(0.0,0.3));
                        }

                        if(no_rotation_bool){
                                // Selecting the verb
                                if(verb < 0.3){
                                        data.move_proposal(0)++;
                                        // Rcpp::Rcout <<  " Grow rotate error" << endl;

                                        grow(all_forest.trees[t],data,partial_residuals,t,
                                                      pow_dist_train,
                                                      pow_dist_test,
                                                      pow_dist_train_test);
                                } else if ((verb>=0.3) & (verb <0.6)) {
                                        data.move_proposal(2)++;
                                        // Rcpp::Rcout <<  " Prune error" << endl;
                                        prune(all_forest.trees[t], data, partial_residuals,t,
                                              pow_dist_train,
                                              pow_dist_test,
                                              pow_dist_train_test);
                                } else {
                                        data.move_proposal(3)++;
                                        change(all_forest.trees[t],data,partial_residuals,t,
                                                        pow_dist_train,
                                                        pow_dist_test,
                                                        pow_dist_train_test);
                                }
                        } else if(only_rotation_bool) {
                                // Selecting the verb
                                if(verb < 0.3){
                                        data.move_proposal(1)++;
                                        // Rcpp::Rcout <<  " Grow rotate error" << endl;

                                        grow_rotation(all_forest.trees[t],data,partial_residuals,t,
                                                      pow_dist_train,
                                                      pow_dist_test,
                                                      pow_dist_train_test);
                                } else if ((verb>=0.3) & (verb <0.6)) {
                                        data.move_proposal(2)++;
                                        // Rcpp::Rcout <<  " Prune error" << endl;
                                        prune(all_forest.trees[t], data, partial_residuals,t,
                                              pow_dist_train,
                                              pow_dist_test,
                                              pow_dist_train_test);
                                } else {
                                        data.move_proposal(4)++;
                                        change_rotation(all_forest.trees[t],data,partial_residuals,t,
                                                        pow_dist_train,
                                                        pow_dist_test,
                                                        pow_dist_train_test);
                                }

                        } else {
                                // Selecting the verb
                                if(verb < 0.15){
                                        // Rcpp::Rcout <<  " Grow error" << endl;
                                        data.move_proposal(0)++;
                                        grow(all_forest.trees[t],data,partial_residuals,t,
                                             pow_dist_train,
                                             pow_dist_test,
                                             pow_dist_train_test);
                                } else if ((verb>=0.15) & (verb< 0.3) ){
                                        data.move_proposal(1)++;
                                        // Rcpp::Rcout <<  " Grow rotate error" << endl;

                                        grow_rotation(all_forest.trees[t],data,partial_residuals,t,
                                                      pow_dist_train,
                                                      pow_dist_test,
                                                      pow_dist_train_test);


                                } else if ((verb>=0.3) & (verb <0.6)) {
                                        data.move_proposal(2)++;
                                        // Rcpp::Rcout <<  " Prune error" << endl;
                                        prune(all_forest.trees[t], data, partial_residuals,t,
                                              pow_dist_train,
                                              pow_dist_test,
                                              pow_dist_train_test);
                                } else if ((verb >= 0.6) & (verb < 0.8)){
                                        data.move_proposal(3)++;
                                        // Rcpp::Rcout <<  " Change error" << endl;
                                        change(all_forest.trees[t], data, partial_residuals,t,
                                               pow_dist_train,
                                               pow_dist_test,
                                               pow_dist_train_test);

                                        // Rcpp::Rcout <<  "Error after change" << endl;
                                } else {
                                        data.move_proposal(4)++;
                                        change_rotation(all_forest.trees[t],data,partial_residuals,t,
                                                        pow_dist_train,
                                                        pow_dist_test,
                                                        pow_dist_train_test);

                                }
                        }

                        // Rcpp::Rcout <<  " Error on phi update" << endl;

                        // Updating phi
                        if(sample_phi){
                                updatePhi(all_forest.trees[t], data, partial_residuals,t,
                                          pow_dist_train,
                                          pow_dist_test,
                                          pow_dist_train_test);
                        }

                        // Getting predictions
                        // Rcpp::Rcout <<  " Error on Get Predictions" << endl;
                        getPredictions(all_forest.trees[t],data,y_hat,prediction_test,t,sample,
                                       y_hat_var,y_hat_test_var,
                                       pow_dist_train,
                                       pow_dist_test,
                                       pow_dist_train_test);

                        // Updating the tree
                        // Rcpp::Rcout <<  "Residuals error 2.0"<< endl;
                        tree_fits_store.col(t) = y_hat;
                        // Rcpp::Rcout <<  "Residuals error 3.0"<< endl;
                        tree_fits_store_test.col(t) = prediction_test;
                        // Rcpp::Rcout <<  "Residuals error 4.0"<< endl;

                        tree_var_store.col(t) = y_hat_var;
                        tree_var_test_store.col(t) = y_hat_test_var;



                }

                // Summing over all trees
                prediction_train_sum = sum(tree_fits_store,1);

                prediction_test_sum = sum(tree_fits_store_test,1);

                var_train_sum = sum(tree_var_store,1);
                var_test_sum = sum(tree_var_test_store,1);

                // Rcpp::Rcout <<  "Error Tau: " << data.tau<< endl;
                updateTau(prediction_train_sum, data);
                // Rcpp::Rcout <<  "New Tau: " << data.tau<< endl;

                // Rcpp::Rcout <<  " All good " << endl;
                tau_post_all(i) = data.tau;

                if(i >= n_burn){
                        // Storing the predictions
                        y_train_hat_post.col(curr) = prediction_train_sum;
                        y_test_hat_post.col(curr) = prediction_test_sum;
                        y_train_var_post.col(curr) = var_train_sum;
                        y_test_var_post.col(curr) = var_test_sum;

                        all_tree_post.slice(curr) = tree_fits_store;
                        all_phi_post.slice(curr) = data.phi_mat;
                        tau_post(curr) = data.tau;
                        curr++;
                }

                if(verbose){
                        pb += 1;
                }

        }

        if(verbose){
                // Initialising PB
                Rcpp::Rcout <<  "[";
                int k = 0;
                // Evaluating progress bar
                for(;k<=pb*width/data.n_mcmc;k++){
                        Rcpp::Rcout <<  "=";
                }

                for(; k < width;k++){
                        Rcpp::Rcout <<  " ";
                }

                Rcpp::Rcout <<  "] " << std::setprecision(5) << 100 << "%\r";
                // std::cout.flush();
                Rcpp::Rcout.flush();

                Rcpp::Rcout <<  std::endl;
        }
        return Rcpp::List::create(y_train_hat_post, // Train prediction [1]
                                  y_test_hat_post, // Test prediction [2]
                                  tau_post, // Tau posterior [3]
                                  all_tree_post, // Each tree posterior [4]
                                  all_phi_post, // All phi posterior [5]
                                  data.move_acceptance, // Acceptance of each verb [6]
                                  data.move_proposal,  // Proposal of each verb [7]
                                  y_train_var_post, // [8]
                                  y_test_var_post, // [9]
                                  tau_post_all // [10]
                                  );
}


//[[Rcpp::export]]
arma::mat mat_init(int n){
        arma::mat A(n,1,arma::fill::ones);
        return A + 4.0;
}


//[[Rcpp::export]]
arma::vec vec_init(int n){
        arma::vec A(n);
        return A+3.0;
}


// Comparing matrix inversions in armadillo
//[[Rcpp::export]]
arma::mat std_inv(arma::mat A, arma::vec diag){

        arma::mat diag_aux = arma::diagmat(diag);
        return arma::inv(A.t()*A+diag_aux);
}

//[[Rcpp::export]]
arma::mat std_pinv(arma::mat A, arma::vec diag){

        arma::mat diag_aux = arma::diagmat(diag);
        return arma::inv_sympd(A.t()*A+diag_aux);
}

//[[Rcpp::export]]
arma::mat faster_simple_std_inv(arma::mat A, arma::vec diag){
        arma::mat diag_aux = arma::diagmat(diag);
        arma::mat L = chol(A.t()*A+diag_aux,"lower");
        return arma::inv(L.t()*L);
}

//[[Rcpp::export]]
double log_test(double a){

        return log(a);
}


//[[Rcpp::export]]
arma::mat faster_std_inv(arma::mat A, arma::vec diag){
        arma::mat ADinvAt = A.t()*arma::diagmat(1.0/diag)*A;
        arma::mat L = arma::chol(ADinvAt + arma::eye(ADinvAt.n_cols,ADinvAt.n_cols),"lower");
        arma::mat invsqrtDA = arma::solve(A.t()/arma::diagmat(arma::sqrt(diag)),L.t());
        arma::mat Ainv = invsqrtDA *invsqrtDA.t()/(ADinvAt + arma::eye(ADinvAt.n_cols,ADinvAt.n_cols));
        return Ainv;
}


//[[Rcpp::export]]
arma::vec rMVN2(const arma::vec& b, const arma::mat& Q)
{
        arma::mat Q_inv = arma::inv(Q);
        arma::mat U = arma::chol(Q_inv, "lower");
        arma::vec z= arma::randn<arma::mat>(Q.n_cols);

        return arma::solve(U.t(), arma::solve(U, z, arma::solve_opts::no_approx), arma::solve_opts::no_approx) + b;
}





//[[Rcpp::export]]
arma::vec rMVNslow(const arma::vec& b, const arma::mat& Q){

        // Rcpp::Rcout <<  "Error sample BETA" << endl;
        arma::vec sample = arma::randn<arma::mat>(Q.n_cols);
        return arma::chol(Q,"lower")*sample + b;

}

//[[Rcpp::export]]
arma::mat matrix_mat(arma::cube array){
        return array.slice(1).t()*array.slice(2);
}









