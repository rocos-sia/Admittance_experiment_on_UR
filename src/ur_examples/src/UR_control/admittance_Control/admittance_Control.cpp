#include <Eigen/Dense>
#include <atomic>
#include <chrono>
#include <kdl/chainiksolvervel_pinv.hpp>
#include <kdl/frames.hpp>
#include <kdl/jntarray.hpp>
#include <plog/Appenders/ColorConsoleAppender.h>
#include <plog/Initializers/RollingFileInitializer.h>
#include <plog/Log.h>
#include <ros/ros.h>
#include <sensor_msgs/JointState.h>
#include <std_msgs/Int64.h>
#include <std_msgs/Int64MultiArray.h>
#include <thread>
#include <trac_ik/trac_ik.hpp>
#include <ur_rtde/rtde_control_interface.h>
#include <ur_rtde/rtde_io_interface.h>
#include <ur_rtde/rtde_receive_interface.h>


using namespace KDL;
using namespace std;
using namespace ur_rtde;
using namespace std::chrono;

class ur_force_control
{
public:
    static const unsigned int joint_num{ 6 };
    //机器人读取变量
private:
    KDL::JntArray joints{ joint_num };
    KDL::JntArray joints_vel{ joint_num };
    KDL::JntArray max_joints_vel{ joint_num };

    std::vector< double > TCP_force{ 0, 0, 0 };
    std::vector< double > force_pos_offset{ 0, 0, 0 };
    std::vector< double > force_vel_offset{ 0, 0, 0 };
    std::vector< double > force_last_vel_offset{ 0, 0, 0 };
    std::vector< double > force_acc_offset{ 0, 0, 0 };
    std::vector< double > force_last_acc_offset{ 0, 0, 0 };

    std::vector< double > TCP_torque{ 0, 0, 0 };
    std::vector< double > torque_pos_offset{ 0, 0, 0 };
    std::vector< double > torque_vel_offset{ 0, 0, 0 };
    std::vector< double > torque_last_vel_offset{ 0, 0, 0 };
    std::vector< double > torque_acc_offset{ 0, 0, 0 };
    std::vector< double > torque_last_acc_offset{ 0, 0, 0 };

    KDL ::Frame reference_frame{ };
    KDL ::Frame current_frame{ };

    std::vector< double > M{ 30, 30, 30, 30, 30, 30 };
    std::vector< double > K{ 150, 150, 150, 100, 100, 100 };
    std::vector< double > B{ 30, 30, 30, 30, 30, 30 };

    double dt{ 0.002 };

    KDL::ChainIkSolverVel_pinv* _ik_vel_ptr;
    ChainFkSolverPos_recursive* _fk_ptr;
    KDL::Twist Cartesian_vel;

    sensor_msgs::JointState command_msg;
    std::shared_ptr< RTDEControlInterface > _ur_control_ptr{ nullptr };

public:
    ur_force_control( KDL::ChainIkSolverVel_pinv* ik_vel_ptr, ChainFkSolverPos_recursive* fk_ptr, std::shared_ptr< RTDEControlInterface > ur_control_ptr )
        : _ik_vel_ptr{ ik_vel_ptr },
          _fk_ptr{ fk_ptr },
          _ur_control_ptr( ur_control_ptr )
    {
        joints( 0 ) = 0 * M_PI / 180;
        joints( 1 ) = ( -145 + 30 ) * M_PI / 180;
        joints( 2 ) = ( -80 - 30 ) * M_PI / 180;
        joints( 3 ) = -45 * M_PI / 180;
        joints( 4 ) = 90 * M_PI / 180;
        joints( 5 ) = 0;

        command_msg.name = ros::V_string{ "shoulder_pan_joint", "shoulder_lift_joint", "elbow_joint", "wrist_1_joint", "wrist_2_joint", "wrist_3_joint" };
        command_msg.position.resize( joint_num );

        for ( int i{ 0 }; i < joint_num; i++ )
            B[ i ] = 2 * 0.2 * sqrt( M[ i ] * K[ i ] );

        _fk_ptr->JntToCart( joints, reference_frame );
        current_frame = reference_frame;

        for ( int i{ 0 }; i < joint_num; i++ )
            max_joints_vel( i ) = 40 * M_PI / 180;
    }

    void calculate_translate( )
    {
        for ( int i{ 0 }; i < 3; i++ )
        {
            force_acc_offset[ i ] = ( TCP_force[ i ] - B[ i ] * force_vel_offset[ i ] - K[ i ] * force_pos_offset[ i ] ) / M[ i ];
            force_vel_offset[ i ] = dt * ( force_acc_offset[ i ] + force_last_acc_offset[ i ] ) / 2 + force_vel_offset[ i ];
            force_pos_offset[ i ] = dt * ( force_vel_offset[ i ] + force_last_vel_offset[ i ] ) / 2 + force_pos_offset[ i ];

            force_last_acc_offset[ i ] = force_acc_offset[ i ];
            force_last_vel_offset[ i ] = force_vel_offset[ i ];
            Cartesian_vel.vel( i )     = force_vel_offset[ i ];

            current_frame.p( i ) = force_pos_offset[ i ];  //更新当前坐标系
        }
    }

    void calculate_rotation( )
    {
        KDL::Vector axis{ };
        axis = ( current_frame.M * reference_frame.M.Inverse( ) ).GetRot( );

        KDL::Vector delta_rot;

        for ( int i{ 0 }; i < 3; i++ )
        {
            torque_pos_offset[ i ] = axis[ i ];
            torque_acc_offset[ i ] = ( TCP_torque[ i ] - B[ i ] * torque_vel_offset[ i ] - K[ i ] * torque_pos_offset[ i ] ) / M[ i ];
            torque_vel_offset[ i ] = dt * ( torque_acc_offset[ i ] + torque_last_acc_offset[ i ] ) / 2 + torque_vel_offset[ i ];

            Cartesian_vel.rot( i ) = torque_vel_offset[ i ];
            delta_rot( i )         = dt * ( torque_vel_offset[ i ] + torque_last_vel_offset[ i ] ) / 2;

            torque_last_acc_offset[ i ] = torque_acc_offset[ i ];
            torque_last_vel_offset[ i ] = torque_vel_offset[ i ];
        }

        current_frame.M = KDL::Rotation::Rot( delta_rot, delta_rot.Norm( ) ) * current_frame.M;
    }

    void calculate( )
    {
        calculate_translate( );
        calculate_rotation( );

        _ik_vel_ptr->CartToJnt( joints, Cartesian_vel, joints_vel );
        KDL::Multiply( joints_vel, dt, joints_vel );
        check_before_move( joints_vel );
        KDL::Add( joints, joints_vel, joints );
    }

    void check_before_move( KDL::JntArray& joints_vel )
    {
        for ( int i{ 0 }; i < joint_num; i++ )
        {
            if ( joints_vel( i ) > max_joints_vel( i ) || TCP_force[ i > 2 ? 2 : i ] > 70 || TCP_torque[ i > 2 ? 2 : i ] > 30 )
            {
                PLOG_ERROR << "joint [" << i << "]  速度过快";
                PLOG_ERROR << "目标速度  = " << joints_vel( i );
                PLOG_ERROR << "允许最大速度  = " << max_joints_vel( i );

                if ( _ur_control_ptr )
                {
                    _ur_control_ptr->servoStop( );
                    _ur_control_ptr->stopScript( );
                }
                exit( 0 );
            }
        }
    }

    void set_force( double force_x, double force_y, double force_z )
    {
        TCP_force[ 0 ] = force_x;
        TCP_force[ 1 ] = force_y;
        TCP_force[ 2 ] = force_z;
        PLOG_DEBUG.printf( "TCP_force  = %f %f %f", TCP_force[ 0 ], TCP_force[ 1 ], TCP_force[ 2 ] );
    }

    void set_torque( double tor_que_x, double tor_que_y, double tor_que_z )
    {
        TCP_torque[ 0 ] = tor_que_x;
        TCP_torque[ 1 ] = tor_que_y;
        TCP_torque[ 2 ] = tor_que_z;
        PLOG_DEBUG.printf( "TCP_torque  = %f %f %f", TCP_torque[ 0 ], TCP_torque[ 1 ], TCP_torque[ 2 ] );
    }

    void set_damp( double value )
    {
        static double damp = 1;
        damp += value;

        for ( int i{ 0 }; i < joint_num; i++ )
            B[ i ] = 2 * damp * sqrt( M[ i ] * K[ i ] );

        PLOG_DEBUG << "damp  = " << damp;
    }

    sensor_msgs::JointState get_joint_msg( )
    {
        for ( int i = 0; i < joint_num; i++ )
        {
            command_msg.position[ i ] = joints.data[ i ];
        }
        return command_msg;
    }

    KDL::JntArray get_joint_vector( )
    {
        return joints;
    }
};

void RunTeleop_damp( const std_msgs::Int64::ConstPtr& msg, ur_force_control* ur_admittance_ptr )
{
    if ( msg->data == 9 )
        ur_admittance_ptr->set_damp( 0.1 );
    else if ( msg->data == 10 )
        ur_admittance_ptr->set_damp( -0.1 );
}

void RunTeleop_force_torque( const std_msgs::Int64MultiArray::ConstPtr& msg, ur_force_control* ur_admittance_ptr )
{
    ur_admittance_ptr->set_force( msg->data[ 0 ], msg->data[ 1 ], msg->data[ 2 ] );
    ur_admittance_ptr->set_torque( msg->data[ 3 ], msg->data[ 4 ], msg->data[ 5 ] );
}

int main( int argc, char* argv[] )
{
    //** 变量初始化 **//
    static plog::ColorConsoleAppender< plog::TxtFormatter > consoleAppender;

    plog::init( plog::debug, &consoleAppender );  // Initialize the logger.
    ros::init( argc, argv, "admittance_Control" );
    ros::NodeHandle nh;
    KDL::Chain UR_chain;
    TRAC_IK::TRAC_IK tracik_solver( "base", "wrist_3_link" );
    std::shared_ptr< RTDEControlInterface > rtde_control_ptr{ nullptr };
    std::shared_ptr< RTDEReceiveInterface > rtde_receive_ptr{ nullptr };

    tracik_solver.getKDLChain( UR_chain );
    ChainFkSolverPos_recursive ur_fk{ UR_chain };
    KDL::ChainIkSolverVel_pinv ik_vel{ UR_chain };

    bool is_remoted_control{ false };
    nh.param( "/admittance_Control/is_remoted_control", is_remoted_control, false );
    if ( is_remoted_control )
    {
        rtde_control_ptr.reset( new RTDEControlInterface{ "192.168.3.101", RTDEControlInterface::FLAG_USE_EXT_UR_CAP } );
        rtde_receive_ptr.reset( new RTDEReceiveInterface{ "192.168.3.101" } );
    }

    ur_force_control ur_admittance{ &ik_vel, &ur_fk, rtde_control_ptr };
    ros::Subscriber sub        = nh.subscribe< std_msgs::Int64 >( "dragging_teleop", 1000, boost::bind( &RunTeleop_damp, _1, &ur_admittance ) );
    ros::Publisher publisher   = nh.advertise< sensor_msgs ::JointState >( "my_joint_state", 1000 );
    ros::Subscriber sub_torque = nh.subscribe< std_msgs::Int64MultiArray >( "force_torque_control", 1000, boost::bind( &RunTeleop_force_torque, _1, &ur_admittance ) );

    ros::AsyncSpinner spinner( 2 );
    spinner.start( );
    ros::Rate hz_500{ 500 };
    //**-------------------------------**//

    if ( is_remoted_control )
    {
        const std::vector< double > joint_q{ 0 * M_PI / 180, ( -145 + 30 ) * M_PI / 180, ( -80 - 30 ) * M_PI / 180, -45 * M_PI / 180, 90 * M_PI / 180, 0 };
        rtde_control_ptr->moveJ( joint_q, 0.3, 0.2, false );  //起始位置初始化}
    }

    KDL::JntArray joint_command_array;
    std::vector< double > joint_command_vector( 6, 0 );
    std::vector< double > tcp_force_torque( 6, 0 );
    std::vector< double > tcp_zero_force_torque( 6, 0 );

    if ( is_remoted_control )
    {
        sleep(3);
        tcp_zero_force_torque = rtde_receive_ptr->getActualTCPForce( );
        PLOG_DEBUG<< tcp_zero_force_torque[0];
        PLOG_DEBUG<< tcp_zero_force_torque[1];
        PLOG_DEBUG<< tcp_zero_force_torque[2];
        sleep(3);


    }

    while ( ros::ok( ) )
    {
        joint_command_array = ur_admittance.get_joint_vector( );
        publisher.publish( ur_admittance.get_joint_msg( ) );

        //** UR 500hz伺服控制 **//
        if ( is_remoted_control )
        {
            double velocity       = 0.5;
            double acceleration   = 0.5;
            double dt             = 1.0 / 500;  // 2ms
            double lookahead_time = 0.1;
            double gain           = 300;

            tcp_force_torque = rtde_receive_ptr->getActualTCPForce( );

            for ( int i = 0; i < ur_force_control::joint_num; i++ )
            {
                joint_command_vector[ i ] = joint_command_array( i );
                tcp_force_torque[ i ]    -= tcp_zero_force_torque[ i ];
            }

            auto t_start = high_resolution_clock::now( );
            rtde_control_ptr->servoJ( joint_command_vector, velocity, acceleration, dt, lookahead_time, gain );

            auto t_stop     = high_resolution_clock::now( );
            auto t_duration = std::chrono::duration< double >( t_stop - t_start );

            if ( t_duration.count( ) < dt )
            {
                std::this_thread::sleep_for( std::chrono::duration< double >( dt - t_duration.count( ) ) );
            }

            ur_admittance.set_force( tcp_force_torque[ 0 ], tcp_force_torque[ 1 ], tcp_force_torque[ 2 ] );
            ur_admittance.set_torque( tcp_force_torque[ 3 ], tcp_force_torque[ 4 ], tcp_force_torque[ 5 ] );
        }

        //**-------------------------------**//

        ur_admittance.calculate( );

        hz_500.sleep( );
    }

    if ( is_remoted_control )
    {
        rtde_control_ptr->servoStop( );
        rtde_control_ptr->stopScript( );
    }

    return 0;
}
