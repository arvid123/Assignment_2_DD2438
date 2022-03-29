using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using System;


namespace UnityStandardAssets.Vehicles.Car
{
    [RequireComponent(typeof(CarController))]
    public class SphereAI4 : MonoBehaviour
    {
        Vector3 GetPathNormal(List<Vector3> slice)
        {
            //public static double FindLinearLeastSquaresFit(List<PointF> points, out double m, out double b)
            // Perform the calculation.
            // Find the values S1, Sx, Sy, Sxx, and Sxy.
            float S1 = slice.Count;
            float Sx = 0;
            float Sz = 0;
            float Sxx = 0;
            float Sxz = 0;
            foreach (Vector3 pt in slice)
            {
                Sx += pt.x;
                Sz += pt.z;
                Sxx += pt.x * pt.x;
                Sxz += pt.x * pt.z;
            }

            // Solve for m and b.
            float m = (Sxz * S1 - Sx * Sz) / (Sxx * S1 - Sx * Sx);
            float b = (Sxz * Sx - Sz * Sxx) / (Sx * Sx - S1 * Sxx);
            Vector3 normal = new Vector3(-1, 0, 1 / m);
            Vector3 forward = slice[1] - slice[0];
            if (Vector3.Cross(forward, normal).y > 0)
            {
                normal = new Vector3(-1, 0, -1 / m);
            }

            Vector2 perp = Vector2.Perpendicular(new Vector2(forward.x, forward.z));

            return new Vector3(perp.x, 0f, perp.y);

        }

        Vector3 GetPathForward(List<Vector3> slice)
        {
            //public static double FindLinearLeastSquaresFit(List<PointF> points, out double m, out double b)
            // Perform the calculation.
            // Find the values S1, Sx, Sy, Sxx, and Sxy.
            float S1 = slice.Count;
            float Sx = 0;
            float Sz = 0;
            float Sxx = 0;
            float Sxz = 0;
            foreach (Vector3 pt in slice)
            {
                Sx += pt.x;
                Sz += pt.z;
                Sxx += pt.x * pt.x;
                Sxz += pt.x * pt.z;
            }

            // Solve for m and b.
            float m = (Sxz * S1 - Sx * Sz) / (Sxx * S1 - Sx * Sx);
            float b = (Sxz * Sx - Sz * Sxx) / (Sx * Sx - S1 * Sxx);

            return (slice[1] - slice[0]).normalized; ;

        }

        Vector3 firstDerivative(Vector3 a, Vector3 b, Vector3 c, float h)
        {
            return (c - a) / (2 * h);
        }

        Vector3 secondDerivative(Vector3 a, Vector3 b, Vector3 c, float h)
        {
            return (a - 2 * b + c) / (h * h);
        }

        float getCurvature(Vector3 a, Vector3 b, Vector3 c, float h)
        {
            var f_p = firstDerivative(a, b, c, h);
            var f_p_p = secondDerivative(a, b, c, h);

            var enumerator = Vector3.Cross(f_p, f_p_p).magnitude;
            var denominator = f_p.sqrMagnitude * f_p.magnitude;

            return enumerator / denominator;
        }

        List<GameObject> FindGameObjectWithLayer(int layer)
        {
            List<GameObject> good_list = new List<GameObject>();
            var allObjects = FindObjectsOfType<GameObject>();
            foreach(GameObject objectToCheck in allObjects){
                if(objectToCheck.layer == layer)
                {
                    good_list.Add(objectToCheck);
                }
            }
            return good_list;
        }

        private CarController m_Car; // the car controller we want to use

        public GameObject terrain_manager_game_object;
        public TerrainManager terrain_manager;

        GameObject leader;
        const int num_followers = 4;
        Vector3[] targets = new Vector3[num_followers];
        Vector3[] target_vels = new Vector3[num_followers];
        Vector3[] friends_pos = new Vector3[num_followers];
        Vector3[] friend_vels = new Vector3[num_followers];
        int[] offsets = new int[]{-1, 1, -2, 2};
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

        float xi_1 = 3f;
        float xi_2 = 5000f; //10000000f;
        float obstacle_radius = 0f;
        float rho_zero = 30f;

        List<GameObject> sphere_list = new List<GameObject>();

        private void Start()
        {
            Time.timeScale = 3f;
            friends = GameObject.FindGameObjectsWithTag("Player");
            // get the car controller
            m_Car = GetComponent<CarController>();
            terrain_manager = terrain_manager_game_object.GetComponent<TerrainManager>();

            friends = GameObject.FindGameObjectsWithTag("Player");
            // get the car controller
            m_Car = GetComponent<CarController>();
            terrain_manager = terrain_manager_game_object.GetComponent<TerrainManager>();

            leader = GameObject.Find("ReplayCar (2)");
            leader_old_pos = leader.transform.position;
            for (int i = 0; i < num_followers; i++)
            {
                targets[i] = leader.transform.position;
            }

            for (int i = 0; i < friends.Length; i++)
            {
                if (friends[i].transform.position.x == transform.position.x)
                {
                    my_id = i;
                }
            }

            float[,] distance_matrix = new float[num_followers, num_followers]; // (i,j) = distance from agent i to point j

            for (int i = 0; i < num_followers; i++)
            {
                for (int j = 0; j < num_followers; j++)
                {
                    distance_matrix[i, j] = (friends[i].transform.position - targets[j]).magnitude;
                }
            }

            for (int l = 0; l < num_followers; l++)
            {
                float min_dist = float.PositiveInfinity;
                int min_i = 0;
                int min_j = 0;
                for (int i = 0; i < num_followers; i++)
                {
                    for (int j = 0; j < num_followers; j++)
                    {
                        if (distance_matrix[i, j] < min_dist)
                        {
                            min_dist = distance_matrix[i, j];
                            min_i = i;
                            min_j = j;
                        }
                    }
                }
                if (min_i == my_id)
                {
                    my_target = min_j;
                    break;
                }
                else
                {
                    for (int k = 0; k < num_followers; k++)
                    {
                        distance_matrix[min_i, k] = float.PositiveInfinity;
                        distance_matrix[k, min_j] = float.PositiveInfinity;
                    }
                }
            }

            Vector3 q_obs = Vector3.positiveInfinity;

            //foreach (GameObject obstacle in GameObject.FindGameObjectsWithTag("cube"))
            foreach (GameObject obstacle in FindGameObjectWithLayer(9))
            {
                Vector3 closest_point = obstacle.GetComponent<Collider>().ClosestPoint(transform.position);
                Vector3 min_distance = closest_point - transform.position;
                if (min_distance.sqrMagnitude < (q_obs - transform.position).sqrMagnitude)
                {
                    q_obs = closest_point;
                }
            }
            
            /*
            // Create potential spheres
            var ter = terrain_manager.myInfo;
            float x_step = (ter.x_high - ter.x_low) / ter.x_N;
            float z_step = (ter.z_high - ter.z_low) / ter.z_N;
            for (int i = 0; i < ter.x_N; i++)
            {
                for (int j = 0; j < ter.z_N; j++)
                {
                    if (ter.traversability[i, j] < 0.5f)
                    {
                        //Create 4 spheres in cell
                        GameObject northwest = GameObject.CreatePrimitive(PrimitiveType.Sphere);
                        GameObject northeast = GameObject.CreatePrimitive(PrimitiveType.Sphere);
                        GameObject southeast = GameObject.CreatePrimitive(PrimitiveType.Sphere);
                        GameObject southwest = GameObject.CreatePrimitive(PrimitiveType.Sphere);

                        northwest.transform.position = new Vector3(ter.get_x_pos(i) - x_step/4, 0.5f, ter.get_z_pos(j) + x_step/4);
                        northeast.transform.position = new Vector3(ter.get_x_pos(i) + x_step / 4, 0.5f, ter.get_z_pos(j) + x_step / 4);
                        southeast.transform.position = new Vector3(ter.get_x_pos(i) + x_step / 4, 0.5f, ter.get_z_pos(j) - x_step / 4);
                        southwest.transform.position = new Vector3(ter.get_x_pos(i) - x_step / 4, 0.5f, ter.get_z_pos(j) - x_step / 4);

                        northwest.layer = 10;
                        northeast.layer = 10;
                        southeast.layer = 10;
                        southwest.layer = 10;

                        northwest.GetComponent<Collider>().enabled = false;
                        northeast.GetComponent<Collider>().enabled = false;
                        southeast.GetComponent<Collider>().enabled = false;
                        southwest.GetComponent<Collider>().enabled = false;

                        List<GameObject> spheres = new List<GameObject>();
                        spheres.Add(northwest);
                        spheres.Add(northeast);
                        spheres.Add(southeast);
                        spheres.Add(southwest);
                        sphere_list.AddRange(spheres);

                        if (northeast.transform.position.x == 185 && northeast.transform.position.z == 305)
                        {
                            Debug.Log("foundit");
                        }

                        foreach (var sphere in spheres)
                        {
                            //foreach (GameObject obstacle in GameObject.FindGameObjectsWithTag("cube"))
                            foreach (GameObject obstacle in FindGameObjectWithLayer(9))
                            {
                                Vector3 closest_point = obstacle.transform.position;
                                Vector3 min_distance = closest_point - sphere.transform.position;
                                if (min_distance.sqrMagnitude < (q_obs - sphere.transform.position).sqrMagnitude)
                                {
                                    q_obs = closest_point;
                                }
                            }
                            Vector3 q_ao = q_obs - sphere.transform.position;
                            float rho_i = Mathf.Abs(q_ao.magnitude - obstacle_radius);
                            float rho_inv = 1f / rho_i;
                            float rho_zero_inv = 1f / rho_zero;
                            float rho_diff = (rho_inv - rho_zero_inv);
                            if (rho_i > rho_zero)
                            {
                                rho_diff = 0;
                            }

                            float u_rep = 0.5f * xi_2 * rho_diff * rho_diff;
                            sphere.transform.localScale = Vector3.one * (u_rep + 1)/1000f;
                            Debug.Log("U_rep " + u_rep);
                        }
                    } else
                    {
                        GameObject obstacle = GameObject.CreatePrimitive(PrimitiveType.Sphere);
                        obstacle.transform.position = new Vector3(ter.get_x_pos(i), 0.5f, ter.get_z_pos(j));
                        obstacle.transform.localScale = obstacle_radius * 2 * Vector3.one;
                        obstacle.GetComponent<Renderer>().material.color = Color.red;

                    }
                }
            }*/

            leader = GameObject.Find("ReplayCar (2)");
            
            switch(transform.position.x) {
                case 239f:
                    my_target = 0;
                    break;
                case 258.3f:
                    my_target = 1;
                    break;
                case 246.8f:
                    my_target = 2;
                    break;
                case 252.6f:
                    my_target = 3;
                    break;
                default:
                    Debug.Log("didn't find target");
                    break;
            }

            my_old_position = transform.position;
            my_old_velocity = Vector3.zero;
            old_ui = Vector3.zero;
        }

        //(2) = 239, 230
        //(3) = 258.3, 230
        //(4) = 246.8, 230
        //(5) = 252.6, 230
        private void FixedUpdate()
        {
            var my_position = transform.position + transform.forward * 2f;
            /*
            foreach (GameObject sphere in sphere_list)
            {
                Vector3 q_att = sphere.transform.position - leader_pos;

                float u_att = Vector3.Dot(q_att, q_att) * 0.5f * xi_1;
                float u_rep = 1000f * sphere.transform.localScale.y - 1;

                sphere.transform.localScale = new Vector3(Mathf.Log(u_att + u_rep*10, 2)/2, (u_rep + 1)/1000f, Mathf.Log(u_att + u_rep*10, 2)/2);
            }*/

            leader_pos = leader.transform.position;
            leader_trajectory.Add(leader_pos);
            leader_back = leader.transform.TransformDirection(Vector3.back);
            leader_left = leader.transform.TransformDirection(Vector3.left);
            int offset = 30;
            for (int i = 0; i < num_followers; i++)
            {
                Vector3 new_pos = new Vector3();
                if (leader_trajectory.Count <= offset * (i + 1))
                {
                    new_pos = leader_trajectory[leader_trajectory.Count - 1];
                    targets[i] = new_pos;
                }
                else
                {
                    new_pos = leader_trajectory[leader_trajectory.Count - offset * (i + 1)];
                    targets[i] = new_pos;

                    List<Vector3> slice = new List<Vector3>();
                    slice = leader_trajectory.GetRange(leader_trajectory.Count - offset * (i + 1), 20);
                    Vector3 path_normal = GetPathNormal(slice);
                    targets[i] += path_normal * offsets[i] * fan_scales[i] * 60f;
                }
                //Vector3 new_pos = leader_pos + leader_back * 10f * (i + 1) + leader_left * offsets[i] * 20f * 0f;

                Debug.DrawLine(my_position, targets[my_target], Color.cyan);

                target_vels[i] = (new_pos - targets[i]) * Time.fixedDeltaTime;

                Debug.DrawLine(new_pos, targets[i], Color.cyan);
            }

            /*if (counter < 40)
            {
                counter++;
                for (int i = 0; i < num_followers; i++)
                {
                    targets[i] = leader.transform.position;
                }
                return;
            }

            leader_pos = leader.transform.position;
            leader_back = leader.transform.TransformDirection(Vector3.back);
            leader_left = leader.transform.TransformDirection(Vector3.left);
            for (int i = 0; i < num_followers; i++) {
                Vector3 new_pos = leader_pos + leader_back * 10f * (i + 1) + leader_left * offsets[i] * 20f;
                Vector3 diff = (new_pos - targets[i]);
                target_vels[i] = (new_pos - targets[i]) / Time.fixedDeltaTime;
                targets[i] = new_pos;
                
            }
            foreach(var target in targets) {
                Debug.DrawLine(leader.transform.position, target, Color.red);
            }*/

            Vector3 q_tar = targets[my_target];
            Vector3 v_tar = target_vels[my_target];
            float theta_tar_signed = Mathf.Deg2Rad * Vector3.SignedAngle(Vector3.right, v_tar, Vector3.down);
            float theta_tar = theta_tar_signed;
            Vector3 q_obs = Vector3.positiveInfinity;

           //foreach (GameObject obstacle in GameObject.FindGameObjectsWithTag("cube"))
           foreach (GameObject obstacle in FindGameObjectWithLayer(9))
           {
                Vector3 closest_point = obstacle.GetComponent<Collider>().ClosestPoint(my_position);
                Vector3 min_distance = closest_point - my_position;    
                if (min_distance.sqrMagnitude < (q_obs - my_position).sqrMagnitude)
                {
                    q_obs = closest_point;
                }
            }
            Vector3 q_ao = q_obs - my_position;

            /*for(int i = 0; i < num_followers; i++)
            {
                Vector3 new_friend_pos = friends[i].transform.position;
                friend_vels[i] = (new_friend_pos - friends_pos[i]) * Time.fixedDeltaTime;
                friends_pos[i] = new_friend_pos;
            }*/

            Vector3 q_at = q_tar - my_position;
            float psi_signed = Mathf.Deg2Rad * Vector3.SignedAngle(Vector3.right, q_at, Vector3.down);
            float psi = psi_signed;
            float rho_i = Mathf.Abs(q_ao.magnitude - obstacle_radius);
            float rho_inv = 1f / rho_i;
            float rho_zero_inv = 1f / rho_zero;
            float rho_diff = (rho_inv - rho_zero_inv);
            if (rho_i > rho_zero)
            {
                rho_diff = 0;
            }
            float mu = (xi_2 * rho_inv * rho_diff) / q_ao.magnitude;
            //float mu = (xi_2 * rho_inv * rho_diff) / q_ao.magnitude <= 0 ? 0.001f : q_ao.magnitude;
            float lambda = (mu * q_ao.magnitude) / (xi_1 * q_at.magnitude);
            float theta_ao_signed = Mathf.Deg2Rad * Vector3.SignedAngle(Vector3.right, q_ao, Vector3.down);
            float theta_ao = theta_ao_signed;
            float psi_bar = Mathf.Atan2((Mathf.Sin(psi) - lambda * Mathf.Sin(theta_ao)) , (Mathf.Cos(psi) - lambda * Mathf.Cos(theta_ao)));
            float term_1 = v_tar.magnitude * Mathf.Cos(theta_tar - psi);
            //Debug.Log("Term 1: " + term_1);
            float term_2 = xi_1 * q_at.magnitude;
            //Debug.Log("Term 2: " + term_2);
            float term_3 = v_tar.sqrMagnitude * Mathf.Sin(theta_tar - psi_bar) * Mathf.Sin(theta_tar - psi_bar);
            //Debug.Log("Term 3: " + term_3 + " v_tar.magnitude " + v_tar.magnitude);
            float norm_v_i = Mathf.Sqrt(Mathf.Pow(v_tar.magnitude * Mathf.Cos(theta_tar - psi) + xi_1 * q_at.magnitude, 2) + v_tar.sqrMagnitude * Mathf.Pow(Mathf.Sin(theta_tar - psi_bar), 2));
            float theta_i = theta_tar;
            if (norm_v_i != 0)
            {
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


            my_old_position = my_position;
            

            float steerAngle = 0f;
            float acceleration = 0f;
            float reversing = 0f;
            float braking = 0f;

            
            Vector3 my_velocity = (my_position - my_old_position) / Time.fixedDeltaTime;
            my_old_position = my_position;
            /*
            Vector3 my_acceleration = (my_velocity - my_old_velocity) / Time.fixedDeltaTime;
            my_old_velocity = my_velocity;
            Vector3 u_i_acceleration = (u_i - old_ui) / Time.fixedDeltaTime;
            old_ui = u_i;


            Vector3 current_position = transform.position;

            //Debug.Log("Distance to goal: " + position_error.magnitude);
            //Debug.DrawLine(transform.position, target_position);

            Vector3 target_velocity = u_i;
            float k_d = 5f;
            float k_p = 1f;
            Vector3 velocity_error = target_velocity - my_velocity;
            Vector3 acceleration_error = u_i_acceleration - my_acceleration;
            Vector3 desired_acceleration = k_p * velocity_error + k_d * acceleration_error;*/


            steerAngle = Vector3.Dot(u_i, transform.right);
            //steerAngle = Vector3.Dot(desired_acceleration, transform.right);
            acceleration = Vector3.Dot(u_i, transform.forward);
            //acceleration = Vector3.Dot(desired_acceleration, transform.forward);

            if (acceleration < 0)
            {
                acceleration = 0.5f; // Fixes reverse direction
            }
            

            RaycastHit hit;
            float maxRange = 50f;
            float diag_margin = 13f;
            //Raycast forward
            if (Physics.Raycast(my_position + transform.up, transform.TransformDirection(Vector3.forward), out hit, maxRange))
            {
                Vector3 closestObstacleInFront = transform.TransformDirection(Vector3.forward) * hit.distance;
                //Debug.DrawRay(transform.position, closestObstacleInFront, Color.yellow);
                //Debug.Log("Did Hit at distance " + hit.distance);
            }
            else
            {
                hit.distance = 51;
            }
            RaycastHit leftHit;
            if (Physics.Raycast(transform.position + transform.up, transform.TransformDirection(Vector3.left), out leftHit, maxRange))
            {
                Vector3 closestObstacleOnLeft = transform.TransformDirection(Vector3.left) * leftHit.distance;
                //Debug.DrawRay(transform.position, closestObstacleOnLeft, Color.red);
                //Debug.Log("Left Hit at distance " + leftHit.distance);
            }
            else
            {
                leftHit.distance = 51;
            }
            //Raycast right
            RaycastHit rightHit;
            if (Physics.Raycast(transform.position + transform.up, transform.TransformDirection(Vector3.right), out rightHit, maxRange))
            {
                Vector3 closestObstacleOnRight = transform.TransformDirection(Vector3.right) * rightHit.distance;
                //Debug.DrawRay(transform.position, closestObstacleOnRight, Color.red);
                //Debug.Log("Right Hit at distance " + rightHit.distance);
            }
            else
            {
                rightHit.distance = 51;
            }

            Vector3 right_diag = transform.TransformDirection(new Vector3(1, 0, 1));
            right_diag.Normalize();
            Ray right_diag_ray = new Ray(my_position, right_diag);
            RaycastHit hitData_right_diag;

            Vector3 left_diag = transform.TransformDirection(new Vector3(-1, 0, 1));
            left_diag.Normalize();
            Ray left_diag_ray = new Ray(my_position, left_diag);
            RaycastHit hitData_left_diag;


            if (Physics.Raycast(right_diag_ray, out hitData_right_diag))
            {
                Vector3 closestObstacleOnRightDiag = right_diag * hitData_right_diag.distance;
                //Debug.DrawRay(transform.position, closestObstacleOnRightDiag, Color.red);
            }
            else
            {
                hitData_right_diag.distance = 51;
            }


            if (Physics.Raycast(left_diag_ray, out hitData_left_diag))
            {
                Vector3 closestObstacleOnLeftDiag = left_diag * hitData_left_diag.distance;
                //Debug.DrawRay(transform.position, closestObstacleOnLeftDiag, Color.red);
            }
            else
            {
                hitData_left_diag.distance = 51;
            }

            if (hitData_left_diag.distance < diag_margin && steerAngle < 0)
            { // Close to left wall
              //acceleration_h = acceleration_h + (desired_distance - hitData_r.distance);

                if (hitData_right_diag.distance < diag_margin && steerAngle > 0)
                { // Close to left wall
                    //acceleration_h = acceleration_h + (desired_distance - hitData_r.distance);
                    Debug.Log("Tight squeeze, don't interrupt turns!");
                    //steerAngle += 0;
                }
                else
                {
                    Debug.Log("Swerve right due to left diag distance " + hitData_left_diag.distance);
                    steerAngle = 0.3f;
                }

            }

            if (hitData_right_diag.distance < diag_margin && steerAngle > 0)
            { // Close to left wall
                //acceleration_h = acceleration_h + (desired_distance - hitData_r.distance);
                Debug.Log("Swerve left due to right diag distance " + hitData_right_diag.distance);
                steerAngle = -0.3f;
            }




            if (backup || hit.distance < 5 || hitData_right_diag.distance < 4 || hitData_left_diag.distance < 4)
            { // Already collided, back up
                Debug.Log("Backing up");
                backup = true;
                if (hit.distance > 8 && hitData_right_diag.distance > 5 && hitData_left_diag.distance > 5)
                {
                    Debug.Log("Backup ended");
                    backup = false;
                }
                reversing = -1f;
                if (Vector3.Dot(my_velocity, transform.forward) < 0)
                {
                    //Moving backwards, invert tires
                    steerAngle = -steerAngle;
                }
                //Debug.Log("Crashed with hit distance " + hit.distance + " and steer angle " + steerAngle);
            }

            float distance_to_target = Vector3.Distance(my_position, targets[my_target]);
            if (distance_to_target < 20f)
            {
                acceleration = 0f;
            }
            if (distance_to_target < 20f)
            {
                braking = 1f;
            }

            //m_Car.Move(steerAngle, Mathf.Clamp(acceleration, 0, 1f), reversing, braking);
            //if (u_i.magnitude > 1f/ Time.fixedDeltaTime)
            //{
            //    u_i = u_i.normalized / Time.fixedDeltaTime;
            //}
            transform.position += u_i * Time.fixedDeltaTime;//*/


            bool fan_in = false;
            bool outer_hit = false;
            bool inner_hit = false;
            //Raycast forward
            List<Vector3> slice_2 = leader_trajectory.GetRange(leader_trajectory.Count - offset * (my_target + 1), 20);
            Vector3 forward = GetPathForward(slice_2);
            if (Physics.Raycast(transform.position + transform.up, transform.forward, out hit, 1000f, LayerMask.GetMask("CubeWalls")))
            {
                Vector3 closestObstacleInFront = transform.forward * hit.distance;
                //Debug.DrawRay(transform.position, closestObstacleInFront, Color.yellow);
                fan_in = hit.distance < 50f && fan_scales[my_target] > 0.01f ? true : fan_in;
                //Debug.Log("Did Hit at distance " + hit.distance);
            }

            // Raycast forward shifted
            Vector3 normal = GetPathNormal(slice_2);
            if (Physics.Raycast(transform.position + transform.up + normal * 5f, transform.forward, out hit, 1000f, LayerMask.GetMask("CubeWalls")))
            {
                Vector3 closestObstacleInFront = transform.forward * hit.distance;
                //Debug.DrawRay(transform.position + normal * 5f, closestObstacleInFront, Color.yellow);
                if (hit.distance < 50f)
                {
                    fan_in = true;
                }

                //Debug.Log("Did Hit at distance " + hit.distance);
            }

            if (Physics.Raycast(transform.position + transform.up - normal * 5f, transform.forward, out hit, 1000f, LayerMask.GetMask("CubeWalls")))
            {
                Vector3 closestObstacleInFront = transform.forward * hit.distance;
                //Debug.DrawRay(transform.position - normal * 5f, closestObstacleInFront, Color.yellow);
                if (hit.distance < 50f)
                {
                    fan_in = true;
                }
                //Debug.Log("Did Hit at distance " + hit.distance);
            }

            if (Physics.Raycast(transform.position + transform.up + normal * 15f, transform.forward, out hit, 1000f, LayerMask.GetMask("CubeWalls")))
            {
                Vector3 closestObstacleInFront = transform.forward * hit.distance;
                //Debug.DrawRay(transform.position + normal * 15f, closestObstacleInFront, Color.yellow);
                if (hit.distance < 50f)
                {
                    outer_hit = true;
                }
                //Debug.Log("Did Hit at distance " + hit.distance);
            }

            if (Physics.Raycast(transform.position + transform.up - normal * 15f, transform.forward, out hit, 1000f, LayerMask.GetMask("CubeWalls")))
            {
                Vector3 closestObstacleInFront = transform.forward * hit.distance;
                //Debug.DrawRay(transform.position - normal * 15f, closestObstacleInFront, Color.yellow);
                if (hit.distance < 50f)
                {
                    inner_hit = true;
                }
                //Debug.Log("Did Hit at distance " + hit.distance);
            }

            Vector3 a = leader_trajectory[leader_trajectory.Count - offset * (my_target + 1) - 10];
            Vector3 b = leader_trajectory[leader_trajectory.Count - offset * (my_target + 1)];
            Vector3 c = leader_trajectory[leader_trajectory.Count - offset * (my_target + 1) + 10];

            var curvature = getCurvature(a, b, c, Time.fixedDeltaTime * 10);

            Debug.Log("Curvature " + curvature);

            /*if (fan_in && fan_scales[my_target] > 0.01f)
            {
                if (!inner_hit)
                {
                    fan_scales[my_target] -= 0.02f;
                }
                else if (!outer_hit)
                {
                    if (fan_scales[my_target] < 2f)
                    {
                        fan_scales[my_target] += 0.02f;
                    }
                }
                else
                {
                    fan_scales[my_target] -= 0.02f;
                }

            }
            else
            {
                if (fan_scales[my_target] < 1f)
                {
                    fan_scales[my_target] += 0.01f;
                }
                else
                {
                    fan_scales[my_target] -= 0.01f;
                }

            }

            if (curvature > 0.01f && fan_scales[my_target] > 0.02f)
            {
                fan_scales[my_target] -= 0.02f;
            }//*/


            //Debug.DrawLine(targets[my_target], targets[my_target] + new Vector3(Mathf.Cos(theta_tar), 0, Mathf.Sin(theta_tar)) * 4f);
            Debug.DrawLine(targets[my_target], targets[my_target] + v_tar, Color.blue);
            Debug.DrawLine(my_position, my_position + new Vector3(Mathf.Cos(theta_i), 0, Mathf.Sin(theta_i)) * norm_v_i);
            //Debug.DrawLine(transform.position, transform.position + u_i, Color.grey);

            Debug.DrawLine(my_position, my_position + new Vector3(Mathf.Cos(theta_ao), 0, Mathf.Sin(theta_ao)) * 10f, Color.green);


        }
    }
}
