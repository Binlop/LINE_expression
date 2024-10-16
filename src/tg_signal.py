#!/usr/bin/env python3

import subprocess
from telegram import Bot
import asyncio
from datetime import datetime
import sys

BOT_TOKEN = '7364323381:AAENOUEpyOUMQkfpaDCyIcW1p54dEvWhTMA'
USER_ID = '762955605'

async def command_signal(command):
    result = run_command(command)
    c = datetime.now()
    current_time = c.strftime('%H:%M:%S')
    if result['success']:
        await send_telegram_message(f'Время: {current_time} \nКоманда *{' '.join(command)}* завершилась')
    else:
        print('Команда завершилась с ошибкой')
        print('Ошибка:', result['error'])
        await send_telegram_message(f'Время: {current_time} \n*Ошибка*: *{''.join(result['error'])}*')


def run_command(command):
    try:
        result = subprocess.run(command, shell=True, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)        
        return {
            'success': True,
            'output': result.stdout,
            'error': None
        }
    except subprocess.CalledProcessError as e:
        return {
            'success': False,
            'output': None,
            'error': f'Команда {e.cmd} завершилась с ошибкой {e.returncode}\n{e.stderr}'
        }

async def send_telegram_message(message):
    bot = Bot(token=BOT_TOKEN)
    await bot.send_message(chat_id=USER_ID, text=message, parse_mode='Markdown')

def main():
    command = sys.argv[1:]
    asyncio.run(command_signal(command))

if __name__ == "__main__":
    main()