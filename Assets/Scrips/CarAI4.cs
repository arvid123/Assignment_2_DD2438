using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using System;


namespace UnityStandardAssets.Vehicles.Car
{
    [RequireComponent(typeof(CarController))]
    public class CarAI4 : MonoBehaviour
    {

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
        TerrainManager terrain_manager;

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
        int my_target;
        Vector3 my_old_position;
        bool backup;
        GameObject[] friends;

        GameObject mySphere;
        Vector3 magic_pos;

        private void Start()
        {
            friends = GameObject.FindGameObjectsWithTag("Player");
            // get the car controller
            m_Car = GetComponent<CarController>();
            terrain_manager = terrain_manager_game_object.GetComponent<TerrainManager>();

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
        }

        //(2) = 239, 230
        //(3) = 258.3, 230
        //(4) = 246.8, 230
        //(5) = 252.6, 230


        private void FixedUpdate()
        {
            leader_pos = leader.transform.position;
            leader_back = leader.transform.TransformDirection(Vector3.back);
            leader_left = leader.transform.TransformDirection(Vector3.left);
            for (int i = 0; i < num_followers; i++) {
                Vector3 new_pos = leader_pos + leader_back * 10f * (i + 1) + leader_left * offsets[i] * 20f;
                target_vels[i] = (new_pos - targets[i]) * Time.fixedDeltaTime;
                targets[i] = new_pos;
            }
            foreach(var target in targets) {
                Debug.DrawLine(leader.transform.position, target, Color.red);
            }

            Vector3 q_tar = targets[my_target];
            Vector3 v_tar = target_vels[my_target];
            float theta_tar_rad_signed = Vector3.SignedAngle(Vector3.right, v_tar, Vector3.up);
            float theta_tar = Mathf.Deg2Rad * theta_tar_rad_signed < 0 ? Mathf.PI + -theta_tar_rad_signed : theta_tar_rad_signed;
            Vector3 q_obs = Vector3.positiveInfinity;

           //foreach (GameObject obstacle in GameObject.FindGameObjectsWithTag("cube"))
           foreach (GameObject obstacle in FindGameObjectWithLayer(9))

                {
                    Vector3 min_distance = obstacle.GetComponent<Collider>().ClosestPoint(transform.position);
                if (min_distance.sqrMagnitude < q_obs.sqrMagnitude)
                {
                    q_obs = min_distance;
                }
            }

            Vector3 q_ao = q_obs - transform.position;

            for(int i = 0; i < num_followers; i++)
            {
                Vector3 new_friend_pos = friends[i].transform.position;
                friend_vels[i] = (new_friend_pos - friends_pos[i]) * Time.fixedDeltaTime;
                friends_pos[i] = new_friend_pos;
            }

            Vector3 q_at = q_tar - transform.position;
            float psi_rad_signed = Vector3.SignedAngle(Vector3.right, q_at, Vector3.up);
            float psi = Mathf.Deg2Rad * psi_rad_signed < 0 ? Mathf.PI + -psi_rad_signed : psi_rad_signed;
            float xi_1 = 2f;
            float xi_2 = 1f;
            float rho_inv = 1f / (q_obs - transform.position).magnitude;
            float rho_zero_inv = 1f / 20f;
            float mu = xi_2 * rho_inv * (rho_inv - rho_zero_inv) / q_ao.magnitude;
            float lambda = (mu * q_ao.magnitude) / (xi_1 * q_at.magnitude);
            float theta_ao_rad_signed = Vector3.SignedAngle(Vector3.right, friend_vels[my_target], Vector3.up);
            float theta_ao = Mathf.Deg2Rad * theta_ao_rad_signed < 0 ? Mathf.PI + theta_ao_rad_signed : theta_ao_rad_signed;
            float psi_bar = Mathf.Atan(Mathf.Abs((Mathf.Sin(psi) - lambda * Mathf.Sin(theta_ao))) / Mathf.Abs((Mathf.Cos(psi) - lambda * Mathf.Cos(theta_ao))));
            float norm_v_i = Mathf.Sqrt(Mathf.Pow(v_tar.magnitude * Mathf.Cos(theta_tar - psi) + xi_1 * q_at.magnitude, 2) + v_tar.sqrMagnitude * Mathf.Sin(theta_tar - psi_bar));
            float theta_i = theta_tar;
            if (norm_v_i != 0)
            {
                theta_i = psi_bar + Mathf.Asin(Mathf.Clamp((v_tar.magnitude * Mathf.Sin(theta_tar - psi_bar) / norm_v_i), -1, 1));
            }

            Vector3 u_i = new Vector3(norm_v_i * Mathf.Cos(theta_i), 0f, norm_v_i * Mathf.Sin(theta_i));

            // PD controller goes here

            float steerAngle = 0f;
            float acceleration = 0f;
            float reversing = 0f;
            float braking = 0f;

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




            //*/


            // keep track of my velocity
            Vector3 my_position = transform.position;
            Vector3 my_velocity = (my_position - my_old_position) / Time.fixedDeltaTime;
            my_old_position = my_position;

            Vector3 current_position = transform.position;
            //Vector3 position_error = target_position - current_position;

            //Debug.Log("Distance to goal: " + position_error.magnitude);
            //Debug.DrawLine(transform.position, target_position);

            Vector3 target_velocity = u_i;
            float k_p = 0f;
            float k_d = 1f;
            Vector3 velocity_error = target_velocity - my_velocity;
            Vector3 desired_acceleration = /*k_p * position_error +*/ k_d * velocity_error;
            desired_acceleration = target_velocity;
            

            steerAngle = Vector3.Dot(desired_acceleration, transform.right);
            acceleration = Vector3.Dot(desired_acceleration, transform.forward);

            /*if (my_velocity.magnitude > 1.5 * hit.distance)
            {
                Debug.Log(my_velocity.magnitude + " speed is bigger than distance " + hit.distance);
                braking = 1;
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

            //Debug.Log("Steer, Accel, Reverse, Braking");
            //Debug.Log(steerAngle + ", " + acceleration + ", " + reversing + ", " + braking);
            m_Car.Move(steerAngle, acceleration, reversing, braking);
        }
    }
}
