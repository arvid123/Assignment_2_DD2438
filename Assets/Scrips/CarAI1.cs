using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using System;
using Assets.Scrips;

namespace UnityStandardAssets.Vehicles.Car
{
    [RequireComponent(typeof(CarController))]
    public class CarAI1 : MonoBehaviour
    {
        private CarController m_Car; // the car controller we want to use

        public GameObject terrain_manager_game_object;
        TerrainManager terrain_manager;

        public GameObject[] friends;
        int k = 3;
        List<Vector3> my_path = new List<Vector3>();
        List<List<Vector3>> paths = new List<List<Vector3>>();
        Vector3 my_old_position;
        bool backup = false;
        float[,] map;
        bool prior_to_root = true;
        Vector3 tree_root_position = Vector3.zero;

        //Function for
        List<TreeNode> ComputeSingleRootCover(Vector3 root_position, List<Edge> edge_list, List<Vector3> vertices, float B) {
            TreeNode root = MST(vertices, edge_list, root_position);
            List<TreeNode> roots = new List<TreeNode>();
            for(int i = 0; i < 3; i++) {
                roots.Add(new TreeNode(root_position, null));
            }

            int a = 0;
            foreach(TreeNode child in root.children) {
                roots[a].children.Add(child);
                a = (a + 1) % 3;
            }


            List<TreeNode> all_roots = new List<TreeNode>();
            List<TreeNode> pseudo_trees = new List<TreeNode>();



            // Step 1 find all medium trees
            // Step 2 if root is light or medium tree we are done
            // Step 3 else get 1 pseudo tree
            // Step 4 if root is light we are done, else go to step 1

            while (true) {
                List<TreeNode> medium_trees = GetMediumTrees(root, B);
                all_roots.AddRange(medium_trees);

                if (medium_trees.Contains(root) || w(root) < B ) {
                    break;
                }

                TreeNode pseudo_tree = GetPseudoTree(root, B);
                if (pseudo_tree != null) {
                    all_roots.Add(pseudo_tree);
                }

                if (w(root) < B) {
                    break;
                }
            }
            if (w(root) < B) {
                all_roots.Add(root);
            }
            return all_roots;

            /*if (all_roots.Count == 3)

            if (medium_trees.Contains(root)) {
                if (medium_trees.Count == 3) {
                    Debug.Log("3 threes with no leftover found");
                    return all_roots;
                }
            }
            else{
                TreeNode pseudo_tree = GetPseudoTree(root, B);
                if (pseudo_tree != null) {
                    pseudo_trees.Add(pseudo_tree);
                    all_roots.Add(pseudo)
                    if(medium_trees.Count + pseudo_trees.Count == 3) {

                    }
                } else {
                    Debug.Log("Root is light tree");
                }
                if (w(root) > B) {
                    return null;
                }
            }

            all_roots.Add(root);
            all_roots.AddRange(medium_trees);
            all_roots.AddRange(pseudo_trees);
            Debug.Log("number of trees: " + all_roots.Count);
            if (all_roots.Count > k + 1) {
                all_roots = null;
            }
            
            return all_roots;*/
        }

        TreeNode GetPseudoTree(TreeNode root, float B) {
            MarkTrees(root, B);
            if (root.weight == Weight.LIGHT) {
                return null;
            }
            TreeNode heavy_light_tree = GetHeavyLightTree(root, B);
            float weight = 0f;
            TreeNode step_dad = new TreeNode(root.position, null);
            foreach(var child in heavy_light_tree.children) {
                weight += w(child);
                step_dad.children.Add(child);
                heavy_light_tree.children.Remove(child);
                if (weight > B) {
                    break;
                }
            }
            return step_dad;
        }

        TreeNode GetHeavyLightTree(TreeNode root, float B) {
            bool has_light_children = true;
            foreach(var child in root.children) {
                if(child.weight != Weight.LIGHT) {
                    has_light_children = false;
                } 
            }
            if (has_light_children) {
                return root;
            } 
            foreach(var child in root.children) {
                TreeNode node = GetHeavyLightTree(child, B);
                if (node != null) {
                    return node;
                }
            }
            return null;
        }

        List<TreeNode> GetMediumTrees(TreeNode root, float B) {
            List<TreeNode> medium_trees = new List<TreeNode>();
            TreeNode medium_tree = FindMediumTree(root, B);
            while (medium_tree != null) {
                medium_trees.Add(medium_tree);
                if (medium_tree == root) {
                    Debug.Log("Perfect tree");
                    break;
                } else {
                    medium_tree.parent.children.Remove(medium_tree);
                }
                medium_tree = FindMediumTree(root, B);
            }
            return medium_trees;
        }

        TreeNode FindMediumTree(TreeNode top_node, float B) {
            float T_e = w(top_node) + 1;
            if (T_e >= B && T_e < 2 * B) {
                return top_node;
            }
            else if (T_e < B) {
                return null;
            }
            TreeNode child_medium_tree = null;
            foreach (TreeNode child in top_node.children) {
                child_medium_tree = FindMediumTree(child, B);
                if (child_medium_tree != null) {
                    return child_medium_tree;
                }
            }
            return null;
        }

        void MarkTrees(TreeNode top_node, float B) {
            float T_e = w(top_node) + 1;
            if (T_e >= B && T_e < 2*B) {
                top_node.weight = Weight.MEDIUM;
            } else if (T_e < B) {
                top_node.weight = Weight.LIGHT;
            } else {
                top_node.weight = Weight.HEAVY;
            }
            if (top_node.children.Count < 1) {
                return;
            }
            foreach(TreeNode child in top_node.children) {
                MarkTrees(child, B);
            }
        }

        float w(TreeNode node) {
            if (node.children.Count < 1) {
                return 0f;
            }

            float sum = 0f;
            foreach(var child in node.children) {
                sum += w(child) + Vector3.Distance(node.position, child.position);
                //sum += w(child) + Mathf.Abs(Vector3.Dot((node.position - child.position), Vector3.right));
            }

            return sum;
        }

        List<TreeNode> KTreeCover(TerrainInfo info, float B) {
            float[,] map = info.traversability;
            List<Edge> edge_list = new List<Edge>();
            List<Vector3> vertices = new List<Vector3>();
            int x_N = info.x_N;
            int z_N = info.z_N;

            // Create edge list
            for(int i = 0; i < x_N-1; i++) {
                for(int j = 0; j < z_N-1; j++) {
                    if (map[i,j] == 0) {
                        Vector3 vertex = new Vector3(info.get_x_pos(i), 0, info.get_z_pos(j));
                        vertices.Add(vertex);
                        if (map[i+1,j] == 0) {
                            Edge right_edge = new Edge(vertex, new Vector3(info.get_x_pos(i + 1), 0, info.get_z_pos(j)));
                            edge_list.Add(right_edge);
                        }
                        if (map[i, j+1] == 0) {
                            Edge up_edge = new Edge(vertex, new Vector3(info.get_x_pos(i), 0, info.get_z_pos(j + 1)));
                            edge_list.Add(up_edge);
                        }
                    }
                }
            }
            // Create root list
            List<Vector3> roots = new List<Vector3>();
            foreach(var friend in friends) {
                roots.Add(new Vector3(info.get_x_pos(info.get_i_index(friend.transform.position.x)), 0f, info.get_z_pos(info.get_j_index(friend.transform.position.z))));
            }

            return ComputeSingleRootCover(roots[0], edge_list, vertices, B);
        }

        List<TreeNode> KTreeCover_new(TerrainInfo info, float B) {
            List<Edge> edge_list = new List<Edge>();
            List<Vector3> vertices = new List<Vector3>();
            int x_N = info.x_N;
            int z_N = info.z_N;
            map = new float[x_N, z_N];
            for(int i = 0; i < x_N; i++) {
                for(int j = 0; j < z_N; j++) {
                    map[i, j] = info.traversability[i, j];
                }
            }
            float x_step = (info.x_high - info.x_low) / x_N;
            float z_step = (info.z_high - info.z_low) / z_N;

            // Create edge list
            for (int i = 0; i < x_N - 1; i++) {
                for (int j = 0; j < z_N - 1; j++) {
                    if (map[i, j] == 0) {
                        bool is_supercell = false;
                        if (map[i + 1, j] == 0 && map[i, j + 1] == 0 && map[i + 1, j + 1] == 0) {
                            is_supercell = true;
                        }

                        //If yes make one vertex with center of 4 as position and mark the squares -1, -2, 
                        //Else make normal cell                                                    -3, -4
                        Vector3 vertex = Vector3.zero;
                        if (is_supercell) {
                            vertex = new Vector3(info.get_x_pos(i) + 0.5f * x_step, 0, info.get_z_pos(j) + 0.5f * z_step);
                            map[i, j] = -1;
                            map[i, j + 1] = -2;
                            map[i + 1, j] = -3;
                            map[i + 1, j + 1] = -4;
                        }
                        else {
                            vertex = new Vector3(info.get_x_pos(i), 0, info.get_z_pos(j));
                        }


                        vertices.Add(vertex);

                        //Adding neighbours: Check (i-1, j), (i, j-1), (i, j+1) 
                        //If their position is 0, add normal neighbour
                        //If their position is <0 add position based on which number it is.
                        if (is_supercell) {
                            switch (map[i, j - 1]) {
                                case 0:
                                    // Normal edge to the left
                                    edge_list.Add(new Edge(vertex, new Vector3(info.get_x_pos(i), 0, info.get_z_pos(j - 1))));
                                    break;
                                case -2:
                                    // Top right of supercell to the left
                                    edge_list.Add(new Edge(vertex, new Vector3(info.get_x_pos(i) + 0.5f * x_step, 0, info.get_z_pos(j - 1) - 0.5f * z_step)));
                                    break;
                                case -4:
                                    // Bottom right of supercell to the left
                                    edge_list.Add(new Edge(vertex, new Vector3(info.get_x_pos(i) - 0.5f * x_step, 0, info.get_z_pos(j - 1) - 0.5f * z_step)));
                                    break;
                                default:
                                    break;
                            }

                            switch (map[i - 1, j]) {
                                case 0:
                                    // Normal edge to the top left
                                    edge_list.Add(new Edge(vertex, new Vector3(info.get_x_pos(i - 1), 0, info.get_z_pos(j))));
                                    break;
                                case -3:
                                    // Bottom left of supercell to the top left
                                    edge_list.Add(new Edge(vertex, new Vector3(info.get_x_pos(i - 1) - 0.5f * x_step, 0, info.get_z_pos(j) + 0.5f * z_step)));
                                    break;
                                case -4:
                                    // Bottom right of supercell to the top left
                                    edge_list.Add(new Edge(vertex, new Vector3(info.get_x_pos(i - 1) - 0.5f * x_step, 0, info.get_z_pos(j) - 0.5f * z_step)));
                                    break;
                                default:
                                    break;
                            }

                            switch (map[i - 1, j + 1]) {
                                case 0:
                                    // Normal edge to the top left
                                    edge_list.Add(new Edge(vertex, new Vector3(info.get_x_pos(i - 1), 0, info.get_z_pos(j + 1))));
                                    break;
                                case -3:
                                    // Bottom left of supercell to the top right
                                    edge_list.Add(new Edge(vertex, new Vector3(info.get_x_pos(i - 1) - 0.5f * x_step, 0, info.get_z_pos(j + 1) + 0.5f * z_step)));
                                    break;
                                default:
                                    break;
                            }

                            switch (map[i, j + 2]) {
                                case -3:
                                    // Bottom left of supercell to the right
                                    edge_list.Add(new Edge(vertex, new Vector3(info.get_x_pos(i) - 0.5f * x_step, 0, info.get_z_pos(j+2) + 0.5f * z_step)));
                                    break;
                                default:
                                    break;
                            }
                        }
                        else {
                            switch (map[i - 1, j]) {
                                case 0:
                                    //Normal cell on top
                                    edge_list.Add(new Edge(vertex, new Vector3(info.get_x_pos(i-1), 0, info.get_z_pos(j))));
                                    break;
                                case -3:
                                    //Bottom left of supercell on top
                                    edge_list.Add(new Edge(vertex, new Vector3(info.get_x_pos(i - 1) - 0.5f * x_step, 0, info.get_z_pos(j) + 0.5f * z_step)));
                                    break;
                                case -4:
                                    //Bottom right of supercell on top
                                    edge_list.Add(new Edge(vertex, new Vector3(info.get_x_pos(i - 1) - 0.5f * x_step, 0, info.get_z_pos(j) - 0.5f * z_step)));
                                    break;
                                default:
                                    break;
                            }

                            switch (map[i, j-1]) {
                                case 0:
                                    //Normal cell to the left
                                    edge_list.Add(new Edge(vertex, new Vector3(info.get_x_pos(i), 0, info.get_z_pos(j - 1))));
                                    break;
                                case -2:
                                    //Top right of supercell on left
                                    edge_list.Add(new Edge(vertex, new Vector3(info.get_x_pos(i) + 0.5f * x_step, 0, info.get_z_pos(j-1) - 0.5f * z_step)));
                                    break;
                                case -4:
                                    //Bottom right of supercell on left
                                    edge_list.Add(new Edge(vertex, new Vector3(info.get_x_pos(i) - 0.5f * x_step, 0, info.get_z_pos(j-1) - 0.5f * z_step)));
                                    break;
                                default:
                                    break;
                            }

                            switch (map[i, j+1]) {
                                case -3:
                                    //Bottom left of supercell on right
                                    edge_list.Add(new Edge(vertex, new Vector3(info.get_x_pos(i) - 0.5f * x_step, 0, info.get_z_pos(j+1) + 0.5f * z_step)));
                                    break;
                                default:
                                    break;
                            }
                        }
                    }
                }
            }
            // Create root list
            List<Vector3> roots = new List<Vector3>();
            foreach (var friend in friends) {
                float closest_distance = 10000f;
                Vector3 current_root = Vector3.zero;
                foreach(Vector3 vertex in vertices) {
                    float dist = Vector3.Distance(vertex, friend.transform.position);
                    if(dist < closest_distance){
                        closest_distance = dist;
                        current_root = vertex;
                    }
                }
                roots.Add(current_root);
            }

            /*foreach(Edge edge in edge_list) {
                Debug.DrawLine(edge.A, edge.B, Color.red, 1000f);
            }
            foreach(Vector3 vertex in vertices) {
                GameObject sphere = GameObject.CreatePrimitive(PrimitiveType.Sphere);
                sphere.transform.position = vertex;
                sphere.transform.localScale = Vector3.one * 5f;
            }*/

            return ComputeSingleRootCover(roots[0], edge_list, vertices, B);
        }

        void DrawTree(TreeNode root, Color line_color) {
            Color color = Color.clear;
            switch (root.weight) {
                case Weight.LIGHT:
                    color = Color.green;
                    break;
                case Weight.MEDIUM:
                    color = Color.yellow;
                    break;
                case Weight.HEAVY:
                    color = Color.red;
                    break;
                default:
                    color = Color.white;
                    break;
            }

            //GameObject sphere = GameObject.CreatePrimitive(PrimitiveType.Sphere);
            //sphere.transform.position = root.position;
            //sphere.transform.localScale = Vector3.one * 10f;
            //sphere.GetComponent<Renderer>().material.color = color;

            foreach (TreeNode child in root.children) {
                Debug.DrawLine(root.position, child.position, line_color, 1000f);
                DrawTree(child, line_color);
            }
        }

        TreeNode MST(List<Vector3> Q, List<Edge> edges, Vector3 root_position) {
            List<Vector3> F = new List<Vector3>();
            List<Edge> EF = new List<Edge>();

            // Create map of "cheapest" edges
            Dictionary<Vector3, Edge> E = new Dictionary<Vector3, Edge>();
            foreach (var vertex in Q) {
                E.Add(vertex, null);
            }
            int total_vertices = Q.Count;
            while(Q.Count > 0){
                Vector3 vertex = Vector3.zero;
                for(int i = 0; i < Q.Count; i++) {
                    vertex = Q[i];
                    if(E[vertex] != null) {
                        break;
                    }
                }
                Q.Remove(vertex);
                F.Add(vertex);
                if (E[vertex] == null) {
                    Debug.Log("Vertex " + vertex.x + ", " + vertex.z + " has no edge");
                }
                else {
                    EF.Add(E[vertex]);
                }

                List<Edge> adjacents = new List<Edge>();
                foreach (Edge edge in edges) {
                    if(edge.A == vertex || edge.B == vertex) {
                        adjacents.Add(edge);
                    }
                }
                foreach (Edge neighbour in adjacents) {
                    if(neighbour.A == vertex) {
                        if (Q.Contains(neighbour.B) && E[neighbour.B] == null){
                            E[neighbour.B] = neighbour;
                        }
                    }
                    else if(neighbour.B == vertex) {
                        if (Q.Contains(neighbour.A) && E[neighbour.A] == null) {
                            E[neighbour.A] = neighbour;
                        }
                    }
                }
            }

            Queue<TreeNode> to_convert = new Queue<TreeNode>();
            TreeNode root = new TreeNode(root_position, null);
            to_convert.Enqueue(root);
            while(to_convert.Count > 0) {
                TreeNode curr_node = to_convert.Dequeue();

                List<Edge> kill_list = new List<Edge>();
                foreach (Edge edge in EF) {
                    if (edge.A == curr_node.position) {
                        TreeNode child = new TreeNode(edge.B, curr_node);
                        curr_node.children.Add(child);
                        to_convert.Enqueue(child);
                        kill_list.Add(edge);
                    }
                    else if(edge.B == curr_node.position) {
                        TreeNode child = new TreeNode(edge.A, curr_node);
                        curr_node.children.Add(child);
                        to_convert.Enqueue(child);
                        kill_list.Add(edge);
                    }
                }
                foreach(Edge victim in kill_list) {
                    EF.Remove(victim);
                }
            }
            return root;
        }

        private void Start()
        {
            Time.timeScale = 6f;
            // get the car controller
            m_Car = GetComponent<CarController>();
            terrain_manager = terrain_manager_game_object.GetComponent<TerrainManager>();
            friends = GameObject.FindGameObjectsWithTag("Player");

            List<TreeNode> k_cover = null;
            float max = 100000f;
            float min = 0f;
            float B = 0f;
            while (max - min > 10f) {
                B = (max - min) / 2f + min;
                List<TreeNode> k_cover_candidate = KTreeCover_new(terrain_manager.myInfo, B);
                int num_trees = k_cover_candidate.Count;
                int allowed_trees = 3;
                if(w(k_cover_candidate[k_cover_candidate.Count - 1]) < B){
                    allowed_trees = 4;
                }
                if(num_trees > allowed_trees) {
                    min = B;
                    Debug.Log("b too low");
                }
                else {
                    max = B;
                }
                Debug.Log("Min: " + min + " Max: " + max);
            }

            float lower_bound = B;

            min = 0f;
            max = 100000f;
            while (max - min > 10f) {
                B = (max - min) / 2f + min;
                List<TreeNode> k_cover_candidate = KTreeCover_new(terrain_manager.myInfo, B);
                int num_trees = k_cover_candidate.Count;
                int allowed_trees = 3;
                if (w(k_cover_candidate[k_cover_candidate.Count - 1]) < B) {
                    allowed_trees = 4;
                }
                if (num_trees < allowed_trees) {
                    max = B;
                    Debug.Log("b too high");
                }
                else {
                    min = B;
                }
                Debug.Log("Min: " + min + " Max: " + max);
            }
            float higher_bound = B;
            Debug.Log("Lower_bound: " + lower_bound + ", higher_bound: " + higher_bound);

            int epochs = 40;
            float step_size = (higher_bound - lower_bound) / epochs;

            List<TreeNode> best_candidate = new List<TreeNode>();
            float min_max = 10000f;
            for (float i = lower_bound; i < higher_bound; i += step_size) {
                List<TreeNode> k_cover_candidate = KTreeCover_new(terrain_manager.myInfo, i);
                float min_weight = 10000f;
                float leftover = 0f;
                foreach (TreeNode tree in k_cover_candidate) {
                    float weight = w(tree);
                    if (weight < B) {
                        leftover = weight;
                    }
                    else if (weight < min_weight) {
                        min_weight = weight;
                    }
                }
                float max_weight = min_weight + leftover;
                if(max_weight < min_max) {
                    min_max = max_weight;
                    best_candidate = k_cover_candidate;
                    Debug.Log("Better candidate found with B = " + i);
                }
            }
            Debug.Log("Num of trees: " + best_candidate.Count);
            foreach(TreeNode tree in best_candidate) {
                Debug.Log(w(tree));
            }//*/


            // TODO: Make car move beside the graph, fix the binary search, and adjust weights to euclidean
            //List<TreeNode> best_candidate = KTreeCover_new(terrain_manager.myInfo, 750f);

            Color[] colors = new Color[] { Color.blue, Color.red, Color.green, Color.yellow, Color.cyan, Color.magenta };
            for (int i = 0; i < best_candidate.Count; i++) {
                //DrawTree(best_candidate[i], colors[i % colors.Length]);
            }

            // Generate path based on each root.
            // Assign self to path (lowest assignment takes leftover)
            for(int i = 0; i < best_candidate.Count; i++) {
                paths.Add(new List<Vector3>());
                GenPath(best_candidate[i], paths[i]);
            }

            my_path = new List<Vector3>();
            Color color = Color.cyan;

            switch(name) {
                case "ArmedCar":
                    PathToRoot(best_candidate[0].position, my_path);
                    my_path.AddRange(paths[0]);
                    tree_root_position = paths[0][0];
                    color = Color.red;
                    break;
                case "ArmedCar (1)":
                    PathToRoot(best_candidate[1].position, my_path);
                    my_path.AddRange(paths[1]);
                    tree_root_position = paths[1][0];
                    color = Color.green;
                    break;
                case "ArmedCar (2)":
                    PathToRoot(best_candidate[2].position, my_path);
                    my_path.AddRange(paths[2]);
                    tree_root_position = paths[2][0];
                    color = Color.blue;
                    break;
                default:
                    break;
            }

            for(int i = 0; i < my_path.Count-2; i++) {
                if (!my_path[i].Equals(my_path[i + 2])) continue;
                if (map[terrain_manager.myInfo.get_i_index(my_path[i + 1].x), terrain_manager.myInfo.get_j_index(my_path[i + 1].z)] == 0) continue;

                Vector3 forward = (my_path[i + 1] - my_path[i]).normalized;
                Vector3 normal = Vector3.Cross(forward, Vector3.up);
                Vector3 up_left = my_path[i + 1] + forward * 5f + normal * 5f;
                Vector3 up_right = my_path[i + 1] + forward * 5f - normal * 5f;
                my_path.Insert(i + 2, up_right);
                my_path.Insert(i + 2, up_left);
            }
            Debug.Log(name + " path count " + my_path.Count);
            Vector3 old_point = my_path[0];
            foreach(Vector3 point in my_path) {
                Debug.DrawLine(old_point, point, color, 10000f);
                old_point = point;
                Debug.Log("drawing");
            }
            
            //*/
        }

        // Add path from starting position to root position.
        void PathToRoot(Vector3 root, List<Vector3> acc){
            List<Vector3> goal = new List<Vector3>();  
            goal.Add(root);
            List<Vector3> start = new List<Vector3>();
            start.Add(transform.position);
            Pathgen pathgen = new Pathgen(terrain_manager, 4f, 0f, "car", goal, start);
            acc.AddRange(pathgen.getBezierPathList(0, 1));
        }

        // Add tree traversal path
        void GenPath(TreeNode root, List<Vector3> acc) {
            acc.Add(root.position);
            
            foreach(TreeNode child in root.children) {
                GenPath(child, acc);
                acc.Add(root.position);
            }
        }

        private void FixedUpdate() {
            // PD controller that follows planned path
            // Reverses if next point is behind (backs out of branches)
            // this is how you access information about the terrain from a simulated laser range finder
            float steerAngle = 0f;
            float acceleration = 0f;
            float reversing = 0f;
            float braking = 0f;
            RaycastHit hit;
            float maxRange = 50f;
            //Raycast forward
            if (Physics.Raycast(transform.position + transform.up, transform.TransformDirection(Vector3.forward), out hit, maxRange)) {
                Vector3 closestObstacleInFront = transform.TransformDirection(Vector3.forward) * hit.distance;
                Debug.DrawRay(transform.position, closestObstacleInFront, Color.yellow);
                //Debug.Log("Did Hit at distance " + hit.distance);
            }
            else {
                hit.distance = 51;
            }

            ///*
            //Raycast left
            RaycastHit leftHit;
            if (Physics.Raycast(transform.position + transform.up, transform.TransformDirection(Vector3.left), out leftHit, maxRange)) {
                Vector3 closestObstacleOnLeft = transform.TransformDirection(Vector3.left) * leftHit.distance;
                Debug.DrawRay(transform.position, closestObstacleOnLeft, Color.red);
                //Debug.Log("Left Hit at distance " + leftHit.distance);
            }
            else {
                leftHit.distance = 51;
            }


            //Raycast right
            RaycastHit rightHit;
            if (Physics.Raycast(transform.position + transform.up, transform.TransformDirection(Vector3.right), out rightHit, maxRange)) {
                Vector3 closestObstacleOnRight = transform.TransformDirection(Vector3.right) * rightHit.distance;
                Debug.DrawRay(transform.position, closestObstacleOnRight, Color.red);
                //Debug.Log("Right Hit at distance " + rightHit.distance);
            }
            else {
                rightHit.distance = 51;
            }

            //45 degree steering

            float diag_margin = 6;
            //Vector3 right_velocity = (my_position - my_old_position) / Time.fixedDeltaTime;
            //Vector3 left_velocity = (my_position - my_old_position) / Time.fixedDeltaTime;

            Vector3 right_diag = transform.TransformDirection(new Vector3(1, 0, 1));
            right_diag.Normalize();
            Ray right_diag_ray = new Ray(transform.position, right_diag);
            RaycastHit hitData_right_diag;

            Vector3 left_diag = transform.TransformDirection(new Vector3(-1, 0, 1));
            left_diag.Normalize();
            Ray left_diag_ray = new Ray(transform.position, left_diag);
            RaycastHit hitData_left_diag;


            if (Physics.Raycast(right_diag_ray, out hitData_right_diag)) {
                Vector3 closestObstacleOnRightDiag = right_diag * hitData_right_diag.distance;
                Debug.DrawRay(transform.position, closestObstacleOnRightDiag, Color.red);
            }
            else {
                hitData_right_diag.distance = 51;
            }


            if (Physics.Raycast(left_diag_ray, out hitData_left_diag)) {
                Vector3 closestObstacleOnLeftDiag = left_diag * hitData_left_diag.distance;
                Debug.DrawRay(transform.position, closestObstacleOnLeftDiag, Color.red);
            }
            else {
                hitData_left_diag.distance = 51;
            }

            // keep track of my velocity
            Vector3 my_position = transform.position;
            Vector3 my_velocity = (my_position - my_old_position) / Time.fixedDeltaTime;
            my_old_position = my_position;

            Vector3 target_position = my_path[0];

            Vector3 current_position = transform.position;
            Vector3 position_error = target_position - current_position;

            //Debug.Log("curr" + current_position);
            //Debug.Log("target" + target_position);
            if(prior_to_root){
                if(my_path[0] == tree_root_position){
                    prior_to_root = false;
                }
            }

            if (position_error.magnitude < 6) {
                //if(Math.Abs(position_error.x) < x_size/2 && Math.Abs(position_error.z) < z_size/2 && my_path.Count() > 1){
                Vector3 previous_point = my_path[0];
                if(my_path.Count > 1) {
                    my_path.RemoveAt(0);
                }
                else { // Route to the start of leftover tree and traverse it
                    tree_root_position = paths[3][0];
                    my_path = new List<Vector3>();
                    PathToRoot(tree_root_position, my_path);
                    my_path.AddRange(paths[3]);
                    prior_to_root = true;
                }
                if (map[terrain_manager.myInfo.get_i_index(my_path[0].x), terrain_manager.myInfo.get_j_index(my_path[0].z)] == 0 || prior_to_root) {
                    target_position = my_path[0];
                } else {
                    Vector3 forward_dir = (my_path[0] - previous_point).normalized;
                    Vector3 normal = Vector3.Cross(forward_dir, Vector3.up);
                    my_path[0] += normal * 12f;
                    target_position = my_path[0];
                }
                //Debug.Log("Removed a point, there are now " + my_path.Count() + " points left!");
            }
            position_error = target_position - current_position;
            //Debug.Log("Distance to goal: " + position_error.magnitude);
            Debug.DrawLine(transform.position, target_position);


            Vector3 target_velocity = Vector3.zero;
            float k_p = 1f;
            float k_d = 1f;
            Vector3 velocity_error = target_velocity - my_velocity;
            Vector3 desired_acceleration = k_p * position_error + k_d * velocity_error;

            acceleration = Vector3.Dot(desired_acceleration, transform.forward);
            if (acceleration < 0) {
                acceleration = 0.2f; // Fixes reverse direction
            }
            steerAngle = Vector3.Dot(desired_acceleration, transform.right);


            if (my_velocity.magnitude > 1.5 * hit.distance) {
                Debug.Log(my_velocity.magnitude + " speed is bigger than distance " + hit.distance);
                braking = 1;
            }

            if (hitData_left_diag.distance < diag_margin && steerAngle < 0) { // Close to left wall
                //acceleration_h = acceleration_h + (desired_distance - hitData_r.distance);

                if (hitData_right_diag.distance < diag_margin && steerAngle > 0) { // Close to left wall
                    //acceleration_h = acceleration_h + (desired_distance - hitData_r.distance);
                    Debug.Log("Tight squeeze, don't interrupt turns!");
                    //steerAngle += 0;
                }
                else {
                    Debug.Log("Swerve right due to left diag distance " + hitData_left_diag.distance);
                    steerAngle = 0.3f;
                }

            }

            if (hitData_right_diag.distance < diag_margin && steerAngle > 0) { // Close to left wall
                //acceleration_h = acceleration_h + (desired_distance - hitData_r.distance);
                Debug.Log("Swerve left due to right diag distance " + hitData_right_diag.distance);
                steerAngle = -0.3f;
            }

            bool is_reversing = Vector3.Angle(transform.forward, my_path[0] - transform.position) > 90;
            bool is_almost_reversing = Vector3.Angle(transform.forward, my_path[0] - transform.position) < 90 && Vector3.Angle(transform.forward, my_path[0] - transform.position) > 80;

            if (is_reversing) {
                reversing = -1;
                if (Vector3.Dot(my_velocity, transform.forward) < 0) {
                    //Moving backwards, invert tires
                    steerAngle = -steerAngle;
                }
            }
            if (is_almost_reversing) {
                acceleration = 1;
            }

            

            if (backup || hit.distance < 5 || hitData_right_diag.distance < 4 || hitData_left_diag.distance < 4) { // Already collided, back up
                Debug.Log("Backing up");
                backup = true;
                if (hit.distance > 8 && hitData_right_diag.distance > 5 && hitData_left_diag.distance > 5) {
                    Debug.Log("Backup ended");
                    backup = false;
                }
                reversing = -1f;
                if (Vector3.Dot(my_velocity, transform.forward) < 0) {
                    //Moving backwards, invert tires
                    steerAngle = -steerAngle;
                }
                //Debug.Log("Crashed with hit distance " + hit.distance + " and steer angle " + steerAngle);
            }


            //Debug.Log("Left diag hit at distance " + hitData_left_diag.distance);
            //Debug.Log("Right diag hit at distance " + hitData_right_diag.distance);


            /*if(leftHit.distance < 7 && steerAngle < 0){
                steerAngle = 0;
            }
            if(rightHit.distance < 7 && steerAngle > 0){
                steerAngle = 0;
            }*/


            /*else if(hit.distance < 20){ // Try to avoid collision
                acceleration = 0.05f * hit.distance;
                //steerAngle = rightHit.distance < leftHit.distance ? -1f : 1f;
                Debug.Log("Getting close, acceleration down to " + acceleration);
            }//*/

            // m_Car.Move(a, b, c, d) is how you control the car
            // a steering [-1, 1] (-1 left, 1 right))
            // b acceleration [0, 1] (0 stand still, 1 full gas ahead)
            // c footbrake [-1, 0], set to -1 makes car reverse at full speed. (non-zero footbrake seems to override acceleration completely)
            // d handbrake [0, 1], brings car to a stop (non-zero handbrake seems to override both acceleration and footbrake, car just stands still)

            Debug.Log("Steer, Accel, Reverse, Braking");
            Debug.Log(steerAngle + ", " + acceleration + ", " + reversing + ", " + braking);
            if (my_path.Count == 1 && position_error.magnitude < 8) {
                Debug.Log("Stopping car.");
                m_Car.Move(0f, 0f, 0f, 1f);
            }
            else {
                m_Car.Move(steerAngle, acceleration, reversing, braking);
            }
        }
    }
}
