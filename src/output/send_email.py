import smtplib, socket
import datetime
import subprocess
from typing import Optional, Union, List
from email.mime.text import MIMEText


def prepare_email(receiver: Union[str, List[str]], sender: str, message: str, subject: Optional[str] = None) -> MIMEText:
    """
    Prepare message to send.
    :param receiver: email address or list of adresses where to send the email
    :param message: message to send
    :param subject: optional subject of the message
    :param sender: address of sender (that should be in From:)
    :return: formatted email message
    """
    # create a Mimetext
    msg = MIMEText(message)
    msg["From"] = sender
    msg["To"] = ", ".join(receiver)
    if subject is None:
        subject = "Snakelines report {date}".format(date=datetime.datetime.now().strftime("%Y-%m-%d %H:%M"))
    msg["Subject"] = subject

    return msg


def send_email_smtp(receiver: Union[str, List[str]], message: str, host: str, port: int, subject: Optional[str] = None,
                    login_name: Optional[str] = None, login_pass: Optional[str] = None, quiet: bool = True) -> bool:
    """
    Sends email using SMTP through host:port gate to receiver with message and subject.
    Additionally, if login name and pass if specified, logins to specified account and sends it from there.
    :param receiver: email address or list of adresses where to send the email
    :param message: message to send
    :param host: host address of email client
    :param port: port of the client
    :param subject: optional subject of the message
    :param login_name: optional login name
    :param login_pass: optional login password
    :param quiet: suppress printouts?
    :return: True if completed successfully
    """
    # setup defaults:
    from_addr = "noreply@noreply.com"
    if login_name is not None:
        from_addr = login_name

    # prepare receiver addresses
    if type(receiver) is str:
        to_addr = [receiver]
    else:
        to_addr = receiver

    # prepare email
    msg = prepare_email(to_addr, from_addr, message, subject).as_string()

    # try to send it
    success = False
    server = None
    try:
        # connect to server
        server = smtplib.SMTP(host=host, port=port)

        # authenticate
        if login_name is not None and login_pass is not None:
            server.starttls()
            server.login(login_name, login_pass)

        # send email
        server.sendmail(from_addr, to_addr, msg)
        success = True

    # catch errors and finalize
    except socket.error as e:
        if not quiet:
            print("Could not connect to {host}:{port:d} - is it listening / up? ({error})".format(host=host, port=port, error=repr(e)))
    except Exception as e:
        if not quiet:
            print("Unknown error ({error})".format(error=repr(e)))
    finally:
        if server is not None:
            server.quit()

    # return success flag
    return success


def send_email_sendmail(receiver: Union[str, List[str]], message: str, subject: Optional[str] = None, sender: Optional[str] = None,
                        command: str = '/usr/sbin/sendmail', params: str = '-t -oi', quiet: bool = True) -> int:
    """
    Sends email using linux command to receiver with message and subject.
    :param receiver: email address or list of adresses where to send the email
    :param message: message to send
    :param subject: optional subject of the message
    :param sender: optional sender address
    :param command: linux command for email sending
    :param params: parameters for the command
    :param quiet: suppress printouts?
    :return: integer returncode of send email command, (0 if success)
    """
    # set defaults
    from_addr = "noreply@noreply.com"
    if sender is not None:
        from_addr = sender

    # prepare receiver addresses
    if type(receiver) is str:
        to_addr = [receiver]
    else:
        to_addr = receiver

    # prepare message
    msg = prepare_email(to_addr, from_addr, message, subject)

    # send message
    p = subprocess.Popen(["{command} {params}".format(command=command, params=params)], stdin=subprocess.PIPE, stderr=subprocess.PIPE, stdout=subprocess.PIPE, shell=True)
    stdout_str, stderr_str = p.communicate(msg.as_bytes())

    # check if successful
    if p.returncode != 0 and not quiet:
        print("Sending of email failed with error (returncode={error}, stderr={stderr}, stdout={stdout}):\n".format(error=p.returncode, stderr=stderr_str, stdout=stdout_str))

    # return returncode
    return p.returncode


# send_email_smtp('mkmarcelgg@gmail.com', 'test', host='smtp.gmail.com', port=587, login_name="snakelines.mailclient@gmail.com", login_pass="hesielko", quiet = False)
# send_email_sendmail('mkmarcelgg@gmail.com', 'test', quiet = False)

