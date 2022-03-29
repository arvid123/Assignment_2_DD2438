using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using System;
using Assets.Scrips;

namespace UnityStandardAssets.Vehicles.Car
{
    [RequireComponent(typeof(CarController))]
    public class Old_CarAI4 : MonoBehaviour
    {

        List<GameObject> FindGameObjectWithLayer(int layer)
        {
            List<GameObject> good_list = new List<GameObject>();
            var allObjects = FindObjectsOfType<GameObject>();
            foreach (GameObject objectToCheck in allObjects)
            {
                if (objectToCheck.layer == layer)
                {
                    good_list.Add(objectToCheck);
                }
            }
            return good_list;
        }

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

            return normal.normalized;

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

        private CarController m_Car; // the car controller we want to use

        public GameObject terrain_manager_game_object;
        TerrainManager terrain_manager;

        GameObject leader;
        const int num_followers = 4;
        Vector3[] targets = new Vector3[num_followers];
        Vector3[] target_vels = new Vector3[num_followers];
        Vector3[] friends_pos = new Vector3[num_followers];
        Vector3[] friend_vels = new Vector3[num_followers];
        int[] offsets = new int[] { -1, 1, -2, 2 };
        Vector3 leader_pos;
        Vector3 leader_back;
        Vector3 leader_left;
        int my_target;
        Vector3 my_old_position;
        Vector3 leader_old_pos;
        List<Vector3> leader_trajectory = new List<Vector3>();
        bool backup;
        GameObject[] friends;
        float[] fan_scales = new float[]{1, 1, 1, 1};
        int my_id;

        private void Start()
        {
            friends = GameObject.FindGameObjectsWithTag("Player");
            // get the car controller
            m_Car = GetComponent<CarController>();
            terrain_manager = terrain_manager_game_object.GetComponent<TerrainManager>();

            leader = GameObject.Find("ReplayCar (2)");
            leader_old_pos = leader.transform.position;
            for(int i = 0; i < num_followers; i++)
            {
                targets[i] = leader.transform.position;
            }

            for(int i = 0; i < friends.Length; i++)
            {
                if(friends[i].transform.position.x == transform.position.x)
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
        }

        //(2) = 239, 230
        //(3) = 258.3, 230
        //(4) = 246.8, 230
        //(5) = 252.6, 230


        private void FixedUpdate()
        {

            leader_pos = leader.transform.position;
            leader_trajectory.Add(leader_pos);
            leader_back = leader.transform.TransformDirection(Vector3.back);
            leader_left = leader.transform.TransformDirection(Vector3.left);
            int offset = 60;
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
                    targets[i] += path_normal * offsets[i] * fan_scales[i] * 20f;
                }
                //Vector3 new_pos = leader_pos + leader_back * 10f * (i + 1) + leader_left * offsets[i] * 20f * 0f;
                

                Debug.DrawLine(transform.position, targets[my_target], Color.cyan);

                

                target_vels[i] = (new_pos - targets[i]) * Time.fixedDeltaTime;

                Debug.DrawLine(new_pos, targets[i], Color.cyan);
            }




            

            // PD controller goes here

            float steerAngle = 0f;
            float acceleration = 0f;
            float reversing = 0f;
            float braking = 0f;

            // keep track of my velocity
            Vector3 target_position = targets[my_target];
            Vector3 my_position = transform.position;
            Vector3 my_velocity = (my_position - my_old_position) / Time.fixedDeltaTime;
            my_old_position = my_position;

            Vector3 current_position = transform.position;

            Vector3 position_error = target_position - current_position;

            //Debug.Log("Distance to goal: " + position_error.magnitude);
            //Debug.DrawLine(transform.position, target_position);

            Vector3 target_velocity = target_vels[my_target];
            float k_p = 2f;
            float k_d = 2f;
            Vector3 velocity_error = target_velocity - my_velocity;
            Vector3 desired_acceleration = k_p * position_error + k_d * velocity_error;


            steerAngle = Vector3.Dot(desired_acceleration, transform.right);
            acceleration = Vector3.Dot(desired_acceleration, transform.forward);



            // this is how you access information about the terrain from a simulated laser range finder
            RaycastHit hit;
            float maxRange = 50f;
            //Raycast forward
            if (Physics.Raycast(transform.position + transform.up, transform.TransformDirection(Vector3.forward), out hit, maxRange))
            {
                Vector3 closestObstacleInFront = transform.TransformDirection(Vector3.forward) * hit.distance;
                Debug.DrawRay(transform.position, closestObstacleInFront, Color.yellow);
                //Debug.Log("Did Hit at distance " + hit.distance);
            }
            else
            {
                hit.distance = 51;
            }


            ///*
            //Raycast left
            RaycastHit leftHit;
            if (Physics.Raycast(transform.position + transform.up, transform.TransformDirection(Vector3.left), out leftHit, maxRange))
            {
                Vector3 closestObstacleOnLeft = transform.TransformDirection(Vector3.left) * leftHit.distance;
                Debug.DrawRay(transform.position, closestObstacleOnLeft, Color.red);
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
                Debug.DrawRay(transform.position, closestObstacleOnRight, Color.red);
                //Debug.Log("Right Hit at distance " + rightHit.distance);
            }
            else
            {
                rightHit.distance = 51;
            }

            //45 degree steering

            float diag_margin = 8;
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


            if (Physics.Raycast(right_diag_ray, out hitData_right_diag))
            {
                Vector3 closestObstacleOnRightDiag = right_diag * hitData_right_diag.distance;
                Debug.DrawRay(transform.position, closestObstacleOnRightDiag, Color.red);
            }
            else
            {
                hitData_right_diag.distance = 51;
            }


            if (Physics.Raycast(left_diag_ray, out hitData_left_diag))
            {
                Vector3 closestObstacleOnLeftDiag = left_diag * hitData_left_diag.distance;
                Debug.DrawRay(transform.position, closestObstacleOnLeftDiag, Color.red);
            }
            else
            {
                hitData_left_diag.distance = 51;
            }


            if (my_velocity.magnitude > 1.5 * hit.distance)
            {
                Debug.Log(my_velocity.magnitude + " speed is bigger than distance " + hit.distance);
                braking = 1;
            }
            if (acceleration < 0)
            {
                acceleration = 0.2f; // Fixes reverse direction
            }

            /*if (hitData_left_diag.distance < diag_margin && steerAngle < 0)
            { // Close to left wall
              //acceleration_h = acceleration_h + (desired_distance - hitData_r.distance);

                if (hitData_right_diag.distance < diag_margin && steerAngle > 0)
                { // Close to left wall
                    //acceleration_h = acceleration_h + (desired_distance - hitData_r.distance);
                    Debug.Log("Tight squeeze, don't interrupt turns!");
                    steerAngle += 0;
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
            }*/

            /*
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
            }*/


            //Debug.Log("Left diag hit at distance " + hitData_left_diag.distance);
            //Debug.Log("Right diag hit at distance " + hitData_right_diag.distance);


            if(leftHit.distance < 7 && steerAngle < 0){
                steerAngle = 0;
            }
            if(rightHit.distance < 7 && steerAngle > 0){
                steerAngle = 0;
            }


            else if(hit.distance < 20){ // Try to avoid collision
                acceleration = 0.05f * hit.distance;
                //steerAngle = rightHit.distance < leftHit.distance ? -1f : 1f;
                Debug.Log("Getting close, acceleration down to " + acceleration);
            }//

            // m_Car.Move(a, b, c, d) is how you control the car
            // a steering [-1, 1] (-1 left, 1 right))
            // b acceleration [0, 1] (0 stand still, 1 full gas ahead)
            // c footbrake [-1, 0], set to -1 makes car reverse at full speed. (non-zero footbrake seems to override acceleration completely)
            // d handbrake [0, 1], brings car to a stop (non-zero handbrake seems to override both acceleration and footbrake, car just stands still)

            //Debug.Log("Steer, Accel, Reverse, Braking");
            //Debug.Log(steerAngle + ", " + acceleration + ", " + reversing + ", " + braking);
            m_Car.Move(steerAngle, acceleration, reversing, braking);
            bool fan_in = false;
            bool outer_hit = false;
            bool inner_hit = false;
            //Raycast forward
            List<Vector3> slice_2 = leader_trajectory.GetRange(leader_trajectory.Count - offset * (my_target + 1), 20);
            Vector3 forward = GetPathForward(slice_2);
            if (Physics.Raycast(transform.position + transform.up, transform.forward, out hit, 1000f, LayerMask.GetMask("CubeWalls")))
            {
                Vector3 closestObstacleInFront = transform.forward * hit.distance;
                Debug.DrawRay(transform.position, closestObstacleInFront, Color.yellow);
                fan_in = hit.distance < 50f && fan_scales[my_target] > 0.01f ? true : fan_in;
                //Debug.Log("Did Hit at distance " + hit.distance);
            }

            // Raycast forward shifted
            Vector3 normal = GetPathNormal(slice_2);
            if (Physics.Raycast(transform.position + transform.up + normal * 5f, transform.forward, out hit, 1000f, LayerMask.GetMask("CubeWalls")))
            {
                Vector3 closestObstacleInFront = transform.forward * hit.distance;
                Debug.DrawRay(transform.position + normal * 5f, closestObstacleInFront, Color.yellow);
                if(hit.distance < 50f)
                {
                    fan_in = true;
                }
                
                //Debug.Log("Did Hit at distance " + hit.distance);
            }

            if (Physics.Raycast(transform.position + transform.up - normal * 5f, transform.forward, out hit, 1000f, LayerMask.GetMask("CubeWalls")))
            {
                Vector3 closestObstacleInFront = transform.forward * hit.distance;
                Debug.DrawRay(transform.position - normal * 5f, closestObstacleInFront, Color.yellow);
                if (hit.distance < 50f)
                {
                    fan_in = true;
                }
                //Debug.Log("Did Hit at distance " + hit.distance);
            }

            if (Physics.Raycast(transform.position + transform.up + normal * 15f, transform.forward, out hit, 1000f, LayerMask.GetMask("CubeWalls")))
            {
                Vector3 closestObstacleInFront = transform.forward * hit.distance;
                Debug.DrawRay(transform.position + normal * 15f, closestObstacleInFront, Color.yellow);
                if (hit.distance < 50f)
                {
                    outer_hit = true;
                }
                //Debug.Log("Did Hit at distance " + hit.distance);
            }

            if (Physics.Raycast(transform.position + transform.up - normal * 15f, transform.forward, out hit, 1000f, LayerMask.GetMask("CubeWalls")))
            {
                Vector3 closestObstacleInFront = transform.forward * hit.distance;
                Debug.DrawRay(transform.position - normal * 15f, closestObstacleInFront, Color.yellow);
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

            if (fan_in && fan_scales[my_target] > 0.01f)
            {
                if (!inner_hit)
                {
                    fan_scales[my_target] -= 0.01f;
                }
                else if (!outer_hit)
                {
                    if(fan_scales[my_target] < 2f)
                    {
                        fan_scales[my_target] += 0.01f;
                    }
                }
                else
                {
                    fan_scales[my_target] -= 0.01f;
                }

            } else
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
            }
        }
    }
}

