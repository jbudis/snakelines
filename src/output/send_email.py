import smtplib, socket
import datetime
import subprocess
from typing import Optional, Union, List


def prepare_email(receiver: Union[str, List[str]], sender: str, message: str, subject: Optional[str] = None) -> str:
    """
    Prepare message to send.
    :param receiver: email address or list of adresses where to send the email
    :param message: message to send
    :param subject: optional subject of the message
    :param sender: address of sender (that should be in From:)
    :return: formatted memail message
    """
    # Add the From: and To: headers at the start!
    email_text = """From: {from_addr}\nTo: {to_addr}\nSubject: {subject}\n\n{message}"""
    return email_text.format(from_addr=sender, to_addr=", ".join(receiver), subject=subject, message=message)


def send_email_smtp(receiver: Union[str, List[str]], message: str, host: str, port: int, subject: Optional[str] = None,
                    login_name: Optional[str] = None, login_pass: Optional[str] = None) -> bool:
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
    :return: True if completed successfully
    """
    # setup defaults:
    from_addr = "noreply@noreply.com"
    if login_name is not None:
        from_addr = login_name

    if type(receiver) is str:
        to_addr = [receiver]
    else:
        to_addr = receiver

    if subject is None:
        subject = "Snakelines report {date}".format(date=datetime.datetime.now().strftime("%Y-%m-%d %H:%M"))

    # prepare email
    msg = prepare_email(to_addr, from_addr, message, subject)

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
        print("Could not connect to {host}:{port:d} - is it listening / up? ({error})".format(host=host, port=port, error=repr(e)))
    except Exception as e:
        print("Unknown error ({error})".format(error=repr(e)))
    finally:
        if server is not None:
            server.quit()

    print(msg, from_addr, to_addr, subject)

    # return success flag
    return success


def send_email_sendmail(receiver: Union[str, List[str]], message: str, subject: Optional[str] = None, sender: Optional[str] = None,
                        command: str = 'sendmail', params: str = '-t -oi') -> int:
    """
    Sends email using linux command to receiver with message and subject.
    :param receiver: email address or list of adresses where to send the email
    :param message: message to send
    :param subject: optional subject of the message
    :param sender: optional sender address
    :param command: linux command for email sending
    :param params: parameters for the command
    :return: integer returncode of send email command, (0 if success)
    """
    # set defaults
    from_addr = "noreply@noreply.com"
    if sender is not None:
        from_addr = sender

    # prepare message
    if type(receiver) is str:
        to_addr = [receiver]
    else:
        to_addr = receiver
    msg = prepare_email(to_addr, from_addr, message, subject)

    # send message
    p = subprocess.Popen("{command} {params}".format(command=command, params=params), stdin=subprocess.PIPE)
    _, stderr = p.communicate(msg)

    # check if successful
    if p.returncode != 0:
        print("Sending of email failed with error:\n" + stderr)

    # return returncode
    return p.returncode


send_email_smtp('mkmarcelgg@gmail.com', 'test', host='smtp.gmail.com', port=587, login_name="snakelines.mailclient@gmail.com", login_pass="hesielko")