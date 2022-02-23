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


        List<List<Vector3>> good_routes = new List<List<Vector3>>(3);

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

            int i = 0;
            for (int j = 0; j < enemy_list.Count; j++)
            {
                partitions[i].Add(friend_list.Count + j);
                i = (i + 1) % 3;
            }

            for(int j = 0; j < 3; j++)
            {

                Pathgen pg = new Pathgen(terrain_manager, 4f, 5f, "car", enemy_list, friend_list);
                var good_route = pg.tspSolver(partitions[j].ToArray());
                good_routes.Add(good_route);
                
                /*var old = good_route[0];
                //Debug.Log(good_route.Count);
                foreach (Vector3 pos in good_route)
                {
                    if(transform.position.x == 210)
                    {
                        Debug.DrawLine(old, pos, Color.blue, 1000f);
                    }
                    else if (transform.position.x == 220)
                    {

                        Debug.DrawLine(old, pos, Color.cyan, 1000f);
                    }
                    else
                    {

                        Debug.DrawLine(old, pos, Color.red, 1000f);
                    }
                    old = pos;
                }*/
            }



            //Debug.Log(String.Format("Goal at {0}, {1}", route[1].x, route[1].z));
            foreach (var good_route in good_routes)
            {
                var old = good_route[0];
                //Debug.Log(good_route.Count);
                foreach (Vector3 pos in good_route)
                {
                    if (good_routes.IndexOf(good_route) == 0)
                    {
                        Debug.DrawLine(old, pos, Color.blue, 1000f);
                    }
                    else if (good_routes.IndexOf(good_route) == 1)
                    {

                        Debug.DrawLine(old, pos, Color.cyan, 1000f);
                    }
                    else
                    {

                        Debug.DrawLine(old, pos, Color.red, 1000f);
                    }
                    //Debug.DrawLine(old, pos, Color.blue, 1000f);
                    old = pos;
                }
            }
        }


        private void FixedUpdate()
        {


            // Execute your path here
            // ...

            Vector3 avg_pos = Vector3.zero;

            foreach (GameObject friend in friends)
            {
                avg_pos += friend.transform.position;
            }
            avg_pos = avg_pos / friends.Length;
            Vector3 direction = (avg_pos - transform.position).normalized;

            bool is_to_the_right = Vector3.Dot(direction, transform.right) > 0f;
            bool is_to_the_front = Vector3.Dot(direction, transform.forward) > 0f;

            float steering = 0f;
            float acceleration = 0;

            if (is_to_the_right && is_to_the_front)
            {
                steering = 1f;
                acceleration = 1f;
            }
            else if (is_to_the_right && !is_to_the_front)
            {
                steering = -1f;
                acceleration = -1f;
            }
            else if (!is_to_the_right && is_to_the_front)
            {
                steering = -1f;
                acceleration = 1f;
            }
            else if (!is_to_the_right && !is_to_the_front)
            {
                steering = 1f;
                acceleration = -1f;
            }

            // this is how you access information about the terrain
            int i = terrain_manager.myInfo.get_i_index(transform.position.x);
            int j = terrain_manager.myInfo.get_j_index(transform.position.z);
            float grid_center_x = terrain_manager.myInfo.get_x_pos(i);
            float grid_center_z = terrain_manager.myInfo.get_z_pos(j);

            Debug.DrawLine(transform.position, new Vector3(grid_center_x, 0f, grid_center_z));


            // this is how you control the car
            //Debug.Log("Steering:" + steering + " Acceleration:" + acceleration);
            m_Car.Move(steering, acceleration, acceleration, 0f);
            //m_Car.Move(0f, -1f, 1f, 0f);


        }
    }
}
