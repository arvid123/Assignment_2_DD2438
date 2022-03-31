using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using System;
using System.Linq;
using Assets.Scrips;

namespace UnityStandardAssets.Vehicles.Car
{
    [RequireComponent(typeof(CarController))]
    public class CarAI5 : MonoBehaviour
    {
        int[,] turret_matrix;

        List<Vector3> GetSafePathToTurret(GameObject[] enemies, TerrainManager terrain_manager, Vector3 start_pos, float acceptance_level) {
            int xN = terrain_manager.myInfo.x_N;
            int zN = terrain_manager.myInfo.z_N;
            float x_low = terrain_manager.myInfo.x_low;
            float x_high = terrain_manager.myInfo.x_high;
            float z_low = terrain_manager.myInfo.z_low;
            float z_high = terrain_manager.myInfo.z_high;
            float x_step = (x_high - x_low) / xN;
            float z_step = (z_high - z_low) / zN;
            int start_i = terrain_manager.myInfo.get_i_index(start_pos.x);
            int start_j = terrain_manager.myInfo.get_j_index(start_pos.z);

            float[,] traversability_copy = new float[xN, zN];
            Array.Copy(terrain_manager.myInfo.traversability, traversability_copy, zN * xN);
            int[,] parent = new int[xN, zN];

            foreach (var enemy in enemies) {
                for (int i = 0; i < xN; i++) {
                    for (int j = 0; j < zN; j++) {
                        if (!(traversability_copy[i, j] == 1)) {
                            float x = terrain_manager.myInfo.get_x_pos(i);
                            float z = terrain_manager.myInfo.get_z_pos(j);
                            float offset_from_center = 0.3f;
                            Vector3 top_right = new Vector3(x + offset_from_center * x_step, 0, z + offset_from_center * z_step);
                            Vector3 top_left = new Vector3(x - offset_from_center * x_step, 0, z + offset_from_center * z_step);
                            Vector3 bot_right = new Vector3(x + offset_from_center * x_step, 0, z - offset_from_center * z_step);
                            Vector3 bot_left = new Vector3(x - offset_from_center * x_step, 0, z - offset_from_center * z_step);

                            int corners_seen = 0;
                            if (!Physics.Linecast(enemy.transform.position, top_right, LayerMask.GetMask("CubeWalls"))) corners_seen++;
                            if (!Physics.Linecast(enemy.transform.position, top_left, LayerMask.GetMask("CubeWalls"))) corners_seen++;
                            if (!Physics.Linecast(enemy.transform.position, bot_right, LayerMask.GetMask("CubeWalls"))) corners_seen++;
                            if (!Physics.Linecast(enemy.transform.position, bot_left, LayerMask.GetMask("CubeWalls"))) corners_seen++;
                            
                            if (corners_seen >= 4) {
                                traversability_copy[i, j] -= 1;
                            }
                        }
                    }
                }
            }
            foreach(var square in danger_squares) {
                int i = terrain_manager.myInfo.get_i_index(square.transform.position.x);
                int j = terrain_manager.myInfo.get_j_index(square.transform.position.z);
                
                switch (traversability_copy[i, j]) {
                    case 0:
                        square.GetComponent<Renderer>().material.color = Color.green;
                        break;
                    case -1:
                        square.GetComponent<Renderer>().material.color = Color.yellow;
                        break;
                    case -2:
                        square.GetComponent<Renderer>().material.color = Color.red;
                        break;
                    default:
                        square.GetComponent<Renderer>().material.color = Color.grey;
                        break;
                }
                //square.GetComponent<Renderer>().material.color = Color.blue;
                //square.SetActive(false);
            }

            Queue<Qnode> Q = new Queue<Qnode>();
            Qnode[,] neighbors = new Qnode[xN, zN];
            Qnode start = new Qnode(start_i, start_j, null);
            neighbors[start_i, start_j] = start;
            Q.Enqueue(start);

            bool found_goal = false;
            Qnode v = null;
            List<Qnode> vs = new List<Qnode>();
            while (Q.Count > 0) {
                v = Q.Dequeue();
                if (traversability_copy[v.i, v.j] == acceptance_level) {
                    found_goal = true;
                    vs.Add(v);
                    continue;
                }
                if (v.j + 1 < zN && traversability_copy[v.i, v.j + 1] < 1 && traversability_copy[v.i, v.j + 1] >= acceptance_level && neighbors[v.i, v.j + 1] == null) {
                    Qnode up = new Qnode(v.i, v.j + 1, v);
                    Q.Enqueue(up);
                    neighbors[up.i, up.j] = up;
                }
                if (v.i + 1 < xN && traversability_copy[v.i + 1, v.j] < 1 && traversability_copy[v.i + 1, v.j] >= acceptance_level  && neighbors[v.i + 1, v.j] == null) {
                    Qnode right = new Qnode(v.i + 1, v.j, v);
                    Q.Enqueue(right);
                    neighbors[right.i, right.j] = right;
                }
                if (v.j - 1 >= 0 && traversability_copy[v.i, v.j - 1] < 1 && traversability_copy[v.i, v.j - 1] >= acceptance_level && neighbors[v.i, v.j - 1] == null) {
                    Qnode down = new Qnode(v.i, v.j - 1, v);
                    Q.Enqueue(down);
                    neighbors[down.i, down.j] = down;
                }
                if (v.i - 1 >= 0 && traversability_copy[v.i - 1, v.j] < 1 && traversability_copy[v.i - 1, v.j] >= acceptance_level && neighbors[v.i - 1, v.j] == null) {
                    Qnode left = new Qnode(v.i - 1, v.j, v);
                    Q.Enqueue(left);
                    neighbors[left.i, left.j] = left;
                }
            }

            float min_distance = float.PositiveInfinity;
            // Find goal with shortest path to nearest turret
            foreach (var node in vs) {
                float distance_to_turret = float.PositiveInfinity;

                foreach (var enemy in enemies) {
                    float x = terrain_manager.myInfo.get_x_pos(node.i);
                    float z = terrain_manager.myInfo.get_z_pos(node.j);
                    Vector3 top_right = new Vector3(x + 0.30f * x_step, 0, z + 0.30f * z_step);
                    Vector3 top_left = new Vector3(x - 0.49f * x_step, 0, z + 0.30f * z_step);
                    Vector3 bot_right = new Vector3(x + 0.49f * x_step, 0, z - 0.49f * z_step);
                    Vector3 bot_left = new Vector3(x - 0.49f * x_step, 0, z - 0.30f * z_step);

                    
                    if (!(Physics.Linecast(enemy.transform.position, top_right, LayerMask.GetMask("CubeWalls")) || Physics.Linecast(enemy.transform.position, top_left, LayerMask.GetMask("CubeWalls")) || Physics.Linecast(enemy.transform.position, bot_left, LayerMask.GetMask("CubeWalls")) || Physics.Linecast(enemy.transform.position, bot_right, LayerMask.GetMask("CubeWalls")))) {
                        distance_to_turret = Vector3.Distance(new Vector3(x, 0, z), enemy.transform.position);
                    }
                }

                if (distance_to_turret < min_distance) {
                    v = node;
                    min_distance = distance_to_turret;
                }
            }

            List<Vector3> path = new List<Vector3>();
            if (found_goal) {
                while (v != null) {
                    path.Add(new Vector3(terrain_manager.myInfo.get_x_pos(v.i), 0, terrain_manager.myInfo.get_z_pos(v.j)));
                    v = v.parent;
                }
            }
            else if (acceptance_level > -5){
                Debug.Log("No solution found");
                return GetSafePathToTurret(enemies, terrain_manager, start_pos, acceptance_level - 1);
            }
            path.Reverse();
            return path;
        }
        List<GameObject> FindGameObjectWithLayer(int layer) {
            List<GameObject> good_list = new List<GameObject>();
            var allObjects = FindObjectsOfType<GameObject>();
            foreach (GameObject objectToCheck in allObjects) {
                if (objectToCheck.layer == layer) {
                    good_list.Add(objectToCheck);
                }
            }
            return good_list;
        }
        Vector3 GetClosestObstaclePoint() {
            Vector3 q_obs = Vector3.positiveInfinity;
            foreach (GameObject obstacle in FindGameObjectWithLayer(9)) {
                Vector3 closest_point = obstacle.transform.position; // GetComponent<Collider>().ClosestPoint(transform.position);
                Vector3 min_distance = closest_point - transform.position;
                if (min_distance.sqrMagnitude < (q_obs - transform.position).sqrMagnitude) {
                    q_obs = closest_point;
                }
            }
            if (!is_leader) {
                Vector3 leader_distance = leader.transform.position - transform.position;
                if (leader_distance.sqrMagnitude < (q_obs - transform.position).sqrMagnitude) {
                    q_obs = leader.transform.position;
                    rho_zero = 5f;
                    obstacle_radius = 3f;
                }
                else {
                    rho_zero = 30f;
                    obstacle_radius = 10.7f;
                }
            }
            return q_obs;
        }
        private CarController m_Car; // the car controller we want to use

        public GameObject terrain_manager_game_object;
        TerrainManager terrain_manager;
        public GameObject game_manager_game_object;
        GameManager game_manager;

        GameObject leader;
        int num_of_turrets = 5;
        const int num_followers = 3;
        Vector3[] targets = new Vector3[num_followers];
        Vector3[] target_vels = new Vector3[num_followers];
        Vector3[] friends_pos = new Vector3[num_followers];
        Vector3[] friend_vels = new Vector3[num_followers];
        int[] offsets = new int[] { -1, 1, -2, 2 };
        Vector3 leader_pos;
        Vector3 leader_back;
        Vector3 leader_left;
        Vector3 leader_old_pos;
        List<Vector3> leader_trajectory = new List<Vector3>();
        float[] fan_scales = new float[] { 1, 1, 1, 1 };
        int my_id;
        int my_target;
        Vector3 my_old_position;
        Vector3 my_old_velocity;
        Vector3 old_ui;
        bool backup;
        GameObject[] friends;
        List<GameObject> repulsive_potential_field = new List<GameObject>();
        List<Vector3> my_path = new List<Vector3>();
        List<float> position_errors = new List<float>();
        List<GameObject> danger_squares = new List<GameObject>();
        int counter;
        float xi_1 = 3f;
        float xi_2 = 5000f; //10000000f;
        float obstacle_radius = 10.7f;
        float rho_zero = 30f;
        bool is_leader = false;
        Vector3[,] old_q_matrix = new Vector3[num_followers, num_followers];

        List<GameObject> sphere_list = new List<GameObject>();

        private void Start() {
            GameObject[] enemies = GameObject.FindGameObjectsWithTag("Enemy");

            List<Vector3> my_pos = new List<Vector3>();
            my_pos.Add(transform.position);


            //Time.timeScale = 2f;
            friends = GameObject.FindGameObjectsWithTag("Player");
            // get the car controller


            leader = GameObject.Find("ArmedCar (4)");
            is_leader = name == "ArmedCar (4)";
            leader_old_pos = leader.transform.position;
            leader_back = -leader.transform.forward;
            leader_left = -leader.transform.right;
            for (int i = 0; i < num_followers; i++) {
                targets[i] = leader.transform.position;
            }

            m_Car = GetComponent<CarController>();
            terrain_manager = terrain_manager_game_object.GetComponent<TerrainManager>();
            game_manager = game_manager_game_object.GetComponent<GameManager>();

            if (is_leader) {
                //Visualize turret vision
                int xN = terrain_manager.myInfo.x_N;
                int zN = terrain_manager.myInfo.z_N;
                float x_low = terrain_manager.myInfo.x_low;
                float x_high = terrain_manager.myInfo.x_high;
                float z_low = terrain_manager.myInfo.z_low;
                float z_high = terrain_manager.myInfo.z_high;
                float x_step = (x_high - x_low) / xN;
                float z_step = (z_high - z_low) / zN;
                for (int i = 0; i < xN; i++) {
                    for (int j = 0; j < zN; j++) {
                        if (terrain_manager.myInfo.traversability[i, j] != 1) {
                            GameObject square = GameObject.CreatePrimitive(PrimitiveType.Cube);
                            square.transform.position = new Vector3(terrain_manager.myInfo.get_x_pos(i), 0.01f, terrain_manager.myInfo.get_z_pos(j));
                            square.transform.localScale = new Vector3(x_step, 0.01f, z_step);
                            danger_squares.Add(square);
                            //Debug.Log("Added a square");
                            Color c = Color.blue;
                            switch (terrain_manager.myInfo.traversability[i, j]) {
                                case 0:
                                    c = Color.green;
                                    break;
                                case -1:
                                    c = Color.yellow;
                                    break;
                                case -2:
                                    c = Color.red;
                                    break;
                                default:
                                    c = Color.grey;
                                    break;
                            }
                            c.a = 0.5f;
                            square.GetComponent<Renderer>().material.color = c;
                            //square.SetActive(false);
                        }
                    }
                }
            
                my_path = GetSafePathToTurret(enemies, terrain_manager, transform.position, -1);
                var old_node = my_path[0];
                foreach (var node in my_path) {
                    Debug.DrawLine(node, old_node, Color.black, 10000f);
                    old_node = node;
                }
            }
            

            for (int i = 0; i < friends.Length; i++) {
                if (friends[i].transform.position.x == transform.position.x) {
                    my_target = i;
                    Debug.Log(name + " " + my_id);
                }
            }
            //for (int i = 0; i < num_followers; i++) {
            //    for (int j = 0; j < num_followers; i++) {
            //        old_q_matrix[i, j] = friends[i].transform.position - friends[j].transform.position;
            //    }
            //}

            

            my_old_position = transform.position;
            my_old_velocity = Vector3.zero;
            old_ui = Vector3.zero;
        }

        //(2) = 239, 230
        //(3) = 258.3, 230
        //(4) = 246.8, 230
        //(5) = 252.6, 230
        private void FixedUpdate() {
            counter += 1;
            Vector3 v_tar = Vector3.zero;
            Vector3 q_tar;

            if (leader == null) {
                leader = GameObject.Find("ArmedCar (3)");
                is_leader = name == "ArmedCar (3)";
                if (leader == null) {
                    leader = GameObject.Find("ArmedCar (2)");
                    is_leader = name == "ArmedCar (2)";
                }

                friends = GameObject.FindGameObjectsWithTag("Player");
            }

            if (is_leader) {
                Debug.Log("New Leader");
                int real_num_of_turrets = GameObject.FindGameObjectsWithTag("Enemy").Count();
                bool turret_has_died = (num_of_turrets != real_num_of_turrets);
                if (turret_has_died) {
                    Debug.Log("Num_of_turrets:" + num_of_turrets + "Game manager says: " + game_manager.number_of_turrets);
                    if (turret_has_died) {
                        num_of_turrets--;
                    }
                    my_path = GetSafePathToTurret(GameObject.FindGameObjectsWithTag("Enemy"), terrain_manager, transform.position, -1);
                    var old_node = my_path[0];
                    foreach (var node in my_path) {
                        Debug.DrawLine(node + Vector3.up, old_node + Vector3.up, Color.cyan, 10000f);
                        old_node = node;
                    }
                }
                

                leader_pos = leader.transform.position;
                leader_back = leader.transform.TransformDirection(Vector3.back);
                leader_left = leader.transform.TransformDirection(Vector3.left);

                foreach (var target in targets) {
                    Debug.DrawLine(leader.transform.position, target, Color.red);
                }

                q_tar = my_path[0];
            }
            else {

                leader_pos = leader.transform.position;
                leader_back = leader.transform.TransformDirection(Vector3.back);
                leader_left = leader.transform.TransformDirection(Vector3.left);

                for (int i = 0; i < num_followers; i++) {
                    Vector3 new_pos = leader_pos - leader_back * 4f + leader_left * (i == 0 ? 4f : -4f);
                    Vector3 diff = (new_pos - targets[i]);
                    target_vels[i] = (new_pos - targets[i]) / Time.fixedDeltaTime;
                    targets[i] = new_pos;

                }

                q_tar = targets[my_target];
                Debug.DrawLine(transform.position, q_tar, Color.cyan);
                v_tar = target_vels[my_target];
            }
            
            //float delta_q_ij = 0;
            //for (int j = 0; j < num_followers; j++) {
            //    if (j != my_target) {
            //        delta_q_ij = ((friends[my_target].transform.position - friends[j].transform.position).magnitude - old_q_matrix[my_target, j].magnitude);
            //    }
            //}

            //for (int i = 0; i < num_followers; i++) {
            //    for (int j = 0; j < num_followers; j++) {
            //        old_q_matrix[i, j] = friends[i].transform.position - friends[j].transform.position;
            //    }
            //}

            //Vector3 delta_q_i = (transform.position - my_old_position);
            //float q_dij = 3f;
            //float r = 2f;
            float theta_tar_signed = Mathf.Deg2Rad * Vector3.SignedAngle(Vector3.right, v_tar, Vector3.down);
            float theta_tar = theta_tar_signed;

            Vector3 q_obs = GetClosestObstaclePoint();
            Vector3 q_ao = q_obs - transform.position;
            Vector3 q_at = q_tar - transform.position;
            float psi_signed = Mathf.Deg2Rad * Vector3.SignedAngle(Vector3.right, q_at, Vector3.down);
            float psi = psi_signed;
            float rho_i = Mathf.Abs(q_ao.magnitude - obstacle_radius);
            float rho_inv = 1f / rho_i;
            float rho_zero_inv = 1f / rho_zero;
            float rho_diff = (rho_inv - rho_zero_inv);
            if (rho_i > rho_zero) {
                rho_diff = 0;
            }
            float mu = (xi_2 * rho_inv * rho_diff) / q_ao.magnitude;
            //float mu = (xi_2 * rho_inv * rho_diff) / q_ao.magnitude <= 0 ? 0.001f : q_ao.magnitude;
            float lambda = (mu * q_ao.magnitude) / (xi_1 * q_at.magnitude);
            float theta_ao_signed = Mathf.Deg2Rad * Vector3.SignedAngle(Vector3.right, q_ao, Vector3.down);
            float theta_ao = theta_ao_signed;
            float psi_bar = Mathf.Atan2((Mathf.Sin(psi) - lambda * Mathf.Sin(theta_ao)), (Mathf.Cos(psi) - lambda * Mathf.Cos(theta_ao)));
            float term_1 = v_tar.magnitude * Mathf.Cos(theta_tar - psi);
            //Debug.Log("Term 1: " + term_1);
            float term_2 = xi_1 * q_at.magnitude;
            //Debug.Log("Term 2: " + term_2);
            float term_3 = v_tar.sqrMagnitude * Mathf.Sin(theta_tar - psi_bar) * Mathf.Sin(theta_tar - psi_bar);
            //Debug.Log("Term 3: " + term_3 + " v_tar.magnitude " + v_tar.magnitude);
            float norm_v_i = Mathf.Sqrt(Mathf.Pow(v_tar.magnitude * Mathf.Cos(theta_tar - psi) + xi_1 * q_at.magnitude, 2) + v_tar.sqrMagnitude * Mathf.Pow(Mathf.Sin(theta_tar - psi_bar), 2));
            float theta_i = theta_tar;
            if (norm_v_i != 0) {
                theta_i = psi_bar + Mathf.Asin(Mathf.Clamp((v_tar.magnitude * Mathf.Sin(theta_tar - psi_bar) / norm_v_i), -1, 1));
            }
            Debug.Log("repulsive potential " + (rho_i <= rho_zero ? 0.5f * xi_2 * Mathf.Abs(rho_diff) * rho_diff : 0));
            //Debug.Log("theta_i " + theta_i * Mathf.Rad2Deg);
            //Debug.Log("psi " + psi);
            //Debug.Log("psi_bar " + psi_bar * Mathf.Rad2Deg);
            //Debug.Log("theta_tar " + theta_tar * Mathf.Rad2Deg);
            //Debug.Log("cos psi " + Mathf.Cos(psi));
            //Debug.Log("sin psi " + Mathf.Sin(psi));
            //Debug.Log("lambda " + lambda);

            

            Vector3 u_i = new Vector3(norm_v_i * Mathf.Cos(theta_i), 0f, norm_v_i * Mathf.Sin(theta_i));
            

            Vector3 my_velocity = (transform.position - my_old_position) / Time.fixedDeltaTime;
            my_old_position = transform.position;
            float steerAngle = 0f;
            float acceleration = 0f;
            float reversing = 0f;
            float braking = 0f;
            int pos_counts = position_errors.Count;
            float k_p = 0.3f;
            float k_d = 1f;
            float k_i = 0.03f;
            Vector3 position_error;
            float velocity_error = 0;
            float integral_error = 0;

            if (is_leader) {
                position_error = my_path[0] - transform.position;
            }
            else {
                position_error = q_tar - transform.position;
            }
            if (pos_counts > 1) {
                velocity_error = (position_error.magnitude - position_errors[pos_counts - 1]) / Time.fixedDeltaTime;
            }
            if (pos_counts > 20) {
                integral_error = position_errors.GetRange(pos_counts - 21, 20).Sum() * Time.fixedDeltaTime;
            }
            acceleration = k_p * position_error.magnitude + k_d * velocity_error + k_i * integral_error;

            Debug.Log("P, D, I: " + k_p * position_error.magnitude + ", " + k_d * velocity_error + ", " + k_i * integral_error);
            position_errors.Add(position_error.magnitude);

            //steerAngle = Vector3.SignedAngle(transform.forward, u_i, Vector3.up) / 45f;
            steerAngle = Vector3.Dot(transform.right, u_i);
            //steerAngle = Vector3.Dot(desired_acceleration, transform.right);


            RaycastHit hit;
            float maxRange = 50f;
            float diag_margin = 13f;
            //Raycast forward
            if (Physics.Raycast(transform.position + transform.up, transform.TransformDirection(Vector3.forward), out hit, maxRange)) {
                Vector3 closestObstacleInFront = transform.TransformDirection(Vector3.forward) * hit.distance;
                Debug.DrawRay(transform.position, closestObstacleInFront, Color.yellow);
                //Debug.Log("Did Hit at distance " + hit.distance);
            }
            else {
                hit.distance = 51;
            }
            RaycastHit leftHit;
            if (Physics.Raycast(transform.position + transform.up, transform.TransformDirection(Vector3.left), out leftHit, maxRange)) {
                Vector3 closestObstacleOnLeft = transform.TransformDirection(Vector3.left) * leftHit.distance;
                //Debug.DrawRay(transform.position, closestObstacleOnLeft, Color.red);
                //Debug.Log("Left Hit at distance " + leftHit.distance);
            }
            else {
                leftHit.distance = 51;
            }
            //Raycast right
            RaycastHit rightHit;
            if (Physics.Raycast(transform.position + transform.up, transform.TransformDirection(Vector3.right), out rightHit, maxRange)) {
                Vector3 closestObstacleOnRight = transform.TransformDirection(Vector3.right) * rightHit.distance;
                //Debug.DrawRay(transform.position, closestObstacleOnRight, Color.red);
                //Debug.Log("Right Hit at distance " + rightHit.distance);
            }
            else {
                rightHit.distance = 51;
            }

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

            if (backup || hit.distance < 5 || hitData_right_diag.distance < 4 || hitData_left_diag.distance < 4) { // Already collided, back up
                Debug.Log("Backing up");
                backup = true;
                if (hit.distance > 7 && hitData_right_diag.distance > 5 && hitData_left_diag.distance > 5) {
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

            if (is_leader) {
                q_tar = my_path[0];
                Debug.DrawLine(transform.position, q_tar, Color.cyan);
                if (Vector3.Distance(transform.position, q_tar) < 10f && my_path.Count > 1) {
                    my_path.RemoveAt(0);
                }
                Debug.Log("Leader velocity: " + my_velocity.magnitude);
                float max_dist = 0f;
                foreach(var friend in friends) {
                    float friend_dist = (friend.transform.position - transform.position).magnitude;
                    max_dist = friend_dist > max_dist ? friend_dist : max_dist;
                }
                if(max_dist > 25f) {
                    braking = 1f;
                    acceleration = 0f;
                }
                else if(max_dist > 10f && my_path.Count < 3){
                    braking = 1f;
                    acceleration = 0f;
                }

                if(my_velocity.magnitude > 100f) {
                    braking = 1f;
                }

                
                m_Car.Move(steerAngle, acceleration, reversing, braking);
            }
            else {
                if (position_error.magnitude < 5f) {
                    braking = 1f;
                    acceleration = 0f;
                    Debug.Log("BREAKING");
                }
                if((leader_pos - transform.position).magnitude < 15f){
                    xi_2 = 500f;
                }
                else {
                    xi_2 = 2000f;
                }
                if(counter > 150) {
                    m_Car.Move(steerAngle, acceleration, reversing, braking);
                }
            }

        //if (u_i.magnitude > 1f/ Time.fixedDeltaTime)
        //{
        //    u_i = u_i.normalized / Time.fixedDeltaTime;
        //}
        //transform.position += u_i * Time.fixedDeltaTime;//*/

        //Debug.DrawLine(targets[my_target], targets[my_target] + new Vector3(Mathf.Cos(theta_tar), 0, Mathf.Sin(theta_tar)) * 4f);
        Debug.DrawLine(targets[my_target], targets[my_target] + v_tar, Color.blue);
        Debug.DrawLine(transform.position, transform.position + new Vector3(Mathf.Cos(theta_i), 0, Mathf.Sin(theta_i)) * norm_v_i);
        //Debug.DrawLine(transform.position, transform.position + u_i, Color.grey);

        Debug.DrawLine(transform.position, transform.position + new Vector3(Mathf.Cos(theta_ao), 0, Mathf.Sin(theta_ao)) * 10f, Color.green);

        }
    }
}
