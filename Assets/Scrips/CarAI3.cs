using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using System;
using Assets.Scrips;

namespace UnityStandardAssets.Vehicles.Car
{
    [RequireComponent(typeof(CarController))]
    public class CarAI3 : MonoBehaviour
    {
        private CarController m_Car; // the car controller we want to use

        public GameObject terrain_manager_game_object;
        TerrainManager terrain_manager;

        static public GameObject[] friends;
        static public GameObject[] enemies;
        List<Vector3> friend_list = new List<Vector3>();
        List<Vector3> enemy_list = new List<Vector3>();
        List<int>[] partitions = new List<int>[3]; 
        List<Vector3> my_path = new List<Vector3>();
        Vector3 my_old_position;
        bool backup = false;

        private void Awake()
        {
            // get the car controller
            m_Car = GetComponent<CarController>();
            terrain_manager = terrain_manager_game_object.GetComponent<TerrainManager>();

            // note that both arrays will have holes when objects are destroyed
            // but for initial planning they should work
            friends = GameObject.FindGameObjectsWithTag("Player");

            for (int i = 0; i < 3; i++)
            {
                partitions[i] = new List<int>();
            }


            for(int i = 0; i < friends.Length; i++)
            {
                var friend_start = friends[i].transform.position;
                friend_list.Add(friend_start);
                partitions[i].Add(i);
            }

            // Plan your path here
            // ...

            
        }
        private void Start()
        {

            enemies = GameObject.FindGameObjectsWithTag("Enemy");
            foreach (var enemy in enemies)
            {
                enemy_list.Add(enemy.transform.position);
            }

            foreach (GameObject obj in enemies)
            {
                //Debug.DrawLine(transform.position, obj.transform.position, Color.black, 10f);
            }

            Pathgen pg1 = new Pathgen(terrain_manager, 4f, 5f, "car", enemy_list, friend_list);
            List<List<Vector3>> good_routes = pg1.vrpSolver();


            //Debug.Log(String.Format("Goal at {0}, {1}", route[1].x, route[1].z));

            foreach (var good_route in good_routes)
            {   

                var old = good_route[0];
                if(old.x == transform.position.x && old.z == transform.position.z){
                    my_path = good_route;
                }

                //Debug.Log(good_route.Count);
                foreach (Vector3 pos in good_route)
                {
                    if (good_routes.IndexOf(good_route) == 0)
                    {
                        UnityEngine.Debug.DrawLine(old, pos, Color.blue, 1000f);
                    }
                    else if (good_routes.IndexOf(good_route) == 1)
                    {

                        UnityEngine.Debug.DrawLine(old, pos, Color.cyan, 1000f);
                    }
                    else
                    {

                        UnityEngine.Debug.DrawLine(old, pos, Color.red, 1000f);
                    }
                    //Debug.DrawLine(old, pos, Color.blue, 1000f);
                    old = pos;
                }
            }
        }


        private void FixedUpdate()
        {


            // this is how you access information about the terrain from the map
            int i = terrain_manager.myInfo.get_i_index(transform.position.x);
            int j = terrain_manager.myInfo.get_j_index(transform.position.z);
            float grid_center_x = terrain_manager.myInfo.get_x_pos(i);
            float grid_center_z = terrain_manager.myInfo.get_z_pos(j);

            float steerAngle = 0f;
            float acceleration = 0f;
            float reversing = 0f;
            float braking = 0f;
            


            Debug.DrawLine(transform.position, new Vector3(grid_center_x, 0f, grid_center_z));

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
            else{
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
            else{
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
            else{
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
            else{
                hitData_right_diag.distance = 51;
            }
            
            
            if (Physics.Raycast(left_diag_ray, out hitData_left_diag))
            {
                Vector3 closestObstacleOnLeftDiag = left_diag * hitData_left_diag.distance;
                Debug.DrawRay(transform.position, closestObstacleOnLeftDiag, Color.red);
            }
            else{
                hitData_left_diag.distance = 51;
            }
            



            //*/

           
            // keep track of my velocity
            Vector3 my_position = transform.position;
            Vector3 my_velocity = (my_position - my_old_position) / Time.fixedDeltaTime;
            my_old_position = my_position;
            
            Vector3 target_position = my_path[0];
            
            Vector3 current_position = transform.position;
            Vector3 position_error = target_position - current_position;

            //Debug.Log("curr" + current_position);
            //Debug.Log("target" + target_position);
            if(position_error.magnitude < 8 && my_path.Count > 1){
            //if(Math.Abs(position_error.x) < x_size/2 && Math.Abs(position_error.z) < z_size/2 && my_path.Count() > 1){
                my_path.RemoveAt(0);
                target_position = my_path[0];
                //Debug.Log("Removed a point, there are now " + my_path.Count() + " points left!");
            }

            //Debug.Log("Distance to goal: " + position_error.magnitude);
            Debug.DrawLine(transform.position, new Vector3(my_path[0].x, 0f, my_path[0].z));


            Vector3 target_velocity = Vector3.zero;
            float k_p = 1f;
            float k_d = 1f;
            Vector3 velocity_error = target_velocity - my_velocity;
            Vector3 desired_acceleration = k_p * position_error + k_d * velocity_error;

            acceleration = Vector3.Dot(desired_acceleration, transform.forward);
            if(acceleration < 0){
                acceleration = 0.2f; // Fixes reverse direction
            }
            steerAngle = Vector3.Dot(desired_acceleration, transform.right);


            if(my_velocity.magnitude > 1.5*hit.distance){
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
                else{
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


            if(backup || hit.distance < 5 || hitData_right_diag.distance < 4 || hitData_left_diag.distance < 4){ // Already collided, back up
                Debug.Log("Backing up");
                backup = true;
                if(hit.distance > 8 && hitData_right_diag.distance > 5 && hitData_left_diag.distance > 5){
                    Debug.Log("Backup ended");
                    backup = false;
                }
                reversing = -1f;
                if(Vector3.Dot(my_velocity, transform.forward) < 0){
                    //Moving backwards, invert tires
                    steerAngle = - steerAngle;
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
            if(my_path.Count == 1 && position_error.magnitude < 8){
                Debug.Log("Stopping car.");
                m_Car.Move(0f, 0f, 0f, 1f);
            }
            else{
                m_Car.Move(steerAngle, acceleration, reversing, braking);
            }
        }
    }
}
